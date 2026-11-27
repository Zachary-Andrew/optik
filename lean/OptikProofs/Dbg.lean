import OptikProofs.Kmer

/-!
# Compacted de Bruijn graphs and partitioned gluing (Definition 7.1, Lemma 7.2, prop:glue)

Both the walker in `include/dbg_unitigs.hpp` and partitioned compaction in `include/dbg_part.hpp` are modelled.

Nodes are oriented k-mers.
`E x X` says the counted data holds the edge from x to X, which overlap in k - 1 bases.
optik's edge masks are symmetric, so `E x X` holds exactly when `E (rc X) (rc x)` does.

Definition 7.1 (`Compacts`) is the walker's rule, all four conditions of it.
x has X as its only successor, X has x as its only predecessor, the two are not one canonical k-mer (no self-loop or hairpin), and neither is a palindrome.
The fourth condition matters only for even k and is missing from the rule's one-line summary in dbg_unitigs.hpp, but the code checks it.

Partitioned compaction (prop:glue) reads one partition's nodes at a time.
Inside a partition it applies `Compacts` directly.
Across partitions each side checks its half of the test and, if that passes, files an open end under the canonical (k+1)-mer of the edge.
`glue_iff_compacts` shows that both halves passing is the same as `Compacts`.
By `rightKey_eq_leftKey` and `key_unique`, the two ends of an edge file under one key, and no other edge shares it.
Together they give `partStep_eq_compacts`, by which partitioned compaction steps across exactly the edges the walker steps across, so both build the same unitigs.
-/

namespace OptikProofs.Dbg

open OptikProofs.Kmer

variable (k : ℕ) (E : List Base → List Base → Prop)

/-- X is x's only successor. -/
def UniqueSucc (x X : List Base) : Prop := E x X ∧ ∀ Y, E x Y → Y = X

/-- x is X's only predecessor. -/
def UniquePred (X x : List Base) : Prop := E x X ∧ ∀ z, E z X → z = x

/-- A k-mer equal to its own reverse complement, possible only for even k. -/
def Palindrome (x : List Base) : Prop := x = rc x

/-- Definition 7.1, the rule the walker applies before stepping from x to X. -/
def Compacts (x X : List Base) : Prop :=
  UniqueSucc E x X ∧ UniquePred E X x ∧ canon X ≠ canon x ∧ ¬Palindrome x ∧ ¬Palindrome X

variable {E}

/-- With symmetric masks, X is x's only successor exactly when rc x's only predecessor is rc X. -/
theorem uniqueSucc_iff_rc (hsym : ∀ x X, E x X ↔ E (rc X) (rc x)) (x X : List Base) :
    UniqueSucc E x X ↔ UniquePred E (rc x) (rc X) := by
  constructor
  · rintro ⟨he, huniq⟩
    refine ⟨(hsym x X).mp he, fun z hz => ?_⟩
    have : E x (rc z) := by simpa using (hsym z (rc x)).mp hz
    rw [← huniq _ this, rc_rc]
  · rintro ⟨he, huniq⟩
    refine ⟨by simpa using (hsym (rc X) (rc x)).mp he, fun Y hY => ?_⟩
    have := huniq (rc Y) ((hsym x Y).mp hY)
    rw [← rc_rc Y, this, rc_rc]

theorem palindrome_rc_iff (x : List Base) : Palindrome (rc x) ↔ Palindrome x := by
  unfold Palindrome
  rw [rc_rc]
  exact eq_comm

/-- The rule is strand-symmetric, since x compacts with X exactly when rc X compacts with rc x.
So walking a unitig backwards on the other strand takes the same steps. -/
theorem compacts_rc (hsym : ∀ x X, E x X ↔ E (rc X) (rc x)) (x X : List Base) :
    Compacts E x X ↔ Compacts E (rc X) (rc x) := by
  unfold Compacts
  rw [uniqueSucc_iff_rc hsym x X, uniqueSucc_iff_rc hsym (rc X) (rc x), rc_rc, rc_rc, canon_rc, canon_rc,
    palindrome_rc_iff, palindrome_rc_iff]
  tauto

/-! ## Gluing across partitions (prop:glue) -/

variable (part : ℕ → ℕ)

/-- A canonical k-mer belongs to one minimizer partition, so both strands of a k-mer share it. -/
def partOf (x : List Base) : ℕ := part (canon x)

/-- As a special case, one canonical k-mer means one partition, so an edge between partitions is never a self-loop or hairpin. -/
theorem canon_ne_of_part_ne {x X : List Base} (h : partOf part x ≠ partOf part X) : canon X ≠ canon x :=
  fun hc => h (by simp [partOf, hc])

/-- The x side of an edge into another partition passed its own half of the test. -/
def OpenRight (x X : List Base) : Prop :=
  UniqueSucc E x X ∧ ¬Palindrome x ∧ partOf part x ≠ partOf part X

/-- For an edge from another partition, the X side passed its half of the test. -/
def OpenLeft (x X : List Base) : Prop :=
  UniquePred E X x ∧ ¬Palindrome X ∧ partOf part x ≠ partOf part X

/-- Across partitions, the two open halves together are exactly Definition 7.1. -/
theorem glue_iff_compacts {x X : List Base} (hcross : partOf part x ≠ partOf part X) :
    OpenRight (E := E) part x X ∧ OpenLeft (E := E) part x X ↔ Compacts E x X := by
  unfold OpenRight OpenLeft Compacts
  have := canon_ne_of_part_ne part hcross
  tauto

/-- A step of partitioned compaction is a compacting edge inside a partition or a glued pair of open ends across two. -/
def PartStep (x X : List Base) : Prop :=
  (partOf part x = partOf part X ∧ Compacts E x X) ∨
    (partOf part x ≠ partOf part X ∧ OpenRight (E := E) part x X ∧ OpenLeft (E := E) part x X)

/-- The walker and partitioned compaction step across exactly the same edges (prop:glue). -/
theorem partStep_eq_compacts : PartStep (E := E) part = Compacts E := by
  funext x X
  apply propext
  unfold PartStep
  by_cases hp : partOf part x = partOf part X
  · simp [hp]
  · rw [glue_iff_compacts part hp]
    simp [hp]

/-! ## Keys of open ends -/

/-- An edge overlaps its k-mers in k - 1 bases. -/
def Overlap (x X : List Base) : Prop := x.length = k ∧ X.length = k ∧ x.tail = X.dropLast

/-- The x side builds its (k+1)-mer from x followed by X's last base. -/
def rightWord (x X : List Base) : List Base := x ++ X.drop (k - 1)

/-- The X side's (k+1)-mer puts x's first base before X. -/
def leftWord (x X : List Base) : List Base := x.take 1 ++ X

variable {k}

theorem rightWord_eq_leftWord {x X : List Base} (h : Overlap k x X) :
    rightWord k x X = leftWord x X := by
  obtain ⟨-, hX, ht⟩ := h
  have hsplit : X = X.dropLast ++ X.drop (k - 1) := by
    rw [List.dropLast_eq_take, hX, List.take_append_drop]
  rw [rightWord, leftWord]
  conv_rhs => rw [hsplit]
  rw [← ht, ← List.append_assoc, ← List.drop_one, List.take_append_drop]

/-- Both ends of an edge file under the same key, the canonical (k+1)-mer, so the glue step can pair them. -/
theorem rightKey_eq_leftKey {x X : List Base} (h : Overlap k x X) :
    canon (rightWord k x X) = canon (leftWord x X) := by
  rw [rightWord_eq_leftWord h]

theorem take_rc (w : List Base) (j : ℕ) : (rc w).take j = rc (w.drop (w.length - j)) := by
  simp [rc, List.take_reverse, List.map_drop]

theorem drop_rc (w : List Base) (j : ℕ) : (rc w).drop j = rc (w.take (w.length - j)) := by
  simp [rc, List.drop_reverse, List.map_take]

/-- An edge's (k+1)-mer gives back its first k-mer as the first k bases and its second k-mer as the last k. -/
theorem rightWord_take_drop (hk : 1 ≤ k) {x X : List Base} (h : Overlap k x X) :
    (rightWord k x X).length = k + 1 ∧ (rightWord k x X).take k = x ∧ (rightWord k x X).drop 1 = X := by
  obtain ⟨hx, hX, ht⟩ := h
  refine ⟨by simp [rightWord, hx, hX]; omega, by simp [rightWord, ← hx], ?_⟩
  rw [rightWord_eq_leftWord ⟨hx, hX, ht⟩, leftWord, List.drop_left' (by simp [hx]; omega)]

/-- Equal keys mean the same bidirected edge, read either way round, so an open end is never glued to an unrelated one. -/
theorem key_unique (hk : 1 ≤ k) {x X y Y : List Base} (h1 : Overlap k x X) (h2 : Overlap k y Y)
    (hkey : canon (rightWord k x X) = canon (rightWord k y Y)) :
    (y = x ∧ Y = X) ∨ (y = rc X ∧ Y = rc x) := by
  obtain ⟨hl1, ht1, hd1⟩ := rightWord_take_drop hk h1
  obtain ⟨hl2, ht2, hd2⟩ := rightWord_take_drop hk h2
  rcases (canon_eq_iff (by rw [hl1, hl2])).mp hkey with heq | heq
  · left
    exact ⟨by rw [← ht2, heq, ht1], by rw [← hd2, heq, hd1]⟩
  · right
    constructor
    · rw [← ht2, heq, take_rc, hl1, show k + 1 - k = 1 by omega, hd1]
    · rw [← hd2, heq, drop_rc, hl1, show k + 1 - 1 = k by omega, ht1]

/-- The claim protocol with min_count > 1.
Each solid k-mer claims its cross-partition edges under their keys, and the glue step drops an edge unless its key is claimed twice.
By `key_unique`, only the edge's own two canonical k-mers can claim its key, so two claims mean both ends are solid and one claim means a non-solid end. -/
theorem claimers_are_endpoints (hk : 1 ≤ k) {x X y Y : List Base} (h1 : Overlap k x X) (h2 : Overlap k y Y)
    (hkey : canon (rightWord k x X) = canon (rightWord k y Y)) :
    (canon y = canon x ∧ canon Y = canon X) ∨ (canon y = canon X ∧ canon Y = canon x) := by
  rcases key_unique hk h1 h2 hkey with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
  · exact Or.inl ⟨rfl, rfl⟩
  · exact Or.inr ⟨canon_rc X, canon_rc x⟩

/-! ## Unitig starts (Lemma 7.2) -/

variable (E)

/-- An oriented k-mer that nothing compacts into, which the walker starts from. -/
def IsStart (z : List Base) : Prop := ¬∃ w, Compacts E w z

/-- A maximal non-cycle unitig, whose consecutive k-mers compact while nothing compacts into its first k-mer and its last compacts with nothing. -/
structure Unitig (p : List (List Base)) : Prop where
  ne_nil : p ≠ []
  steps : ∀ i (h : i + 1 < p.length), Compacts E p[i] p[i + 1]
  first : IsStart E (p.head ne_nil)
  last : ¬∃ w, Compacts E (p.getLast ne_nil) w

variable {E}

/-- Lemma 7.2.
Among the oriented forms of a non-cycle unitig's k-mers, the walker's starts are exactly its first k-mer and the reverse complement of its last.
Each unitig therefore has two starts (one, for a palindromic singleton), and emitting it only from the smaller one gives it exactly once. -/
theorem unitig_starts (hsym : ∀ x X, E x X ↔ E (rc X) (rc x)) {p : List (List Base)} (hu : Unitig E p)
    {z : List Base} (hz : z ∈ p ∨ z ∈ p.map rc) :
    IsStart E z ↔ z = p.head hu.ne_nil ∨ z = rc (p.getLast hu.ne_nil) := by
  have hlen : 0 < p.length := List.length_pos_iff.mpr hu.ne_nil
  constructor
  · intro hstart
    rcases hz with hz | hz
    · obtain ⟨i, hi, rfl⟩ := List.getElem_of_mem hz
      left
      rcases Nat.eq_zero_or_pos i with rfl | hpos
      · simp [List.head_eq_getElem]
      · exact absurd ⟨p[i - 1], by simpa [show i - 1 + 1 = i by omega] using hu.steps (i - 1) (by omega)⟩ hstart
    · obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hz
      obtain ⟨i, hi, rfl⟩ := List.getElem_of_mem hy
      right
      by_cases hlast : i + 1 < p.length
      · have := (compacts_rc hsym _ _).mp (hu.steps i hlast)
        exact absurd ⟨_, this⟩ hstart
      · simp [List.getLast_eq_getElem, show i = p.length - 1 by omega]
  · rintro (rfl | rfl)
    · exact hu.first
    · rintro ⟨w, hw⟩
      rw [compacts_rc hsym, rc_rc] at hw
      exact hu.last ⟨_, hw⟩

/-- Walking a unitig's other start gives the same k-mers on the other strand, in reverse order. -/
theorem unitig_rc (hsym : ∀ x X, E x X ↔ E (rc X) (rc x)) {p : List (List Base)} (hu : Unitig E p) :
    Unitig E (p.reverse.map rc) := by
  have hne : p.reverse.map rc ≠ [] := by simpa using hu.ne_nil
  refine ⟨hne, fun i h => ?_, ?_, ?_⟩
  · simp only [List.length_map, List.length_reverse] at h
    simp only [List.getElem_map, List.getElem_reverse]
    have := (compacts_rc hsym _ _).mp (hu.steps (p.length - 1 - (i + 1)) (by omega))
    simpa [show p.length - 1 - (i + 1) + 1 = p.length - 1 - i by omega] using this
  · rintro ⟨w, hw⟩
    rw [List.head_map, List.head_reverse, compacts_rc hsym, rc_rc] at hw
    exact hu.last ⟨_, hw⟩
  · rintro ⟨w, hw⟩
    rw [List.getLast_map, List.getLast_reverse, compacts_rc hsym] at hw
    exact hu.first ⟨_, by simpa using hw⟩

end OptikProofs.Dbg
