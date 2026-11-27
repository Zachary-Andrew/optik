import OptikProofs.Spss
import OptikProofs.Dbg

/-!
# Partitioned counting and edge masks

Models pass 1 and pass 2 of `count_into_graph` in `src/kmer_count.cpp`.

Pass 1 gives every k-mer the minimum order value of its m-mers and cuts each run of k-mers into records, each a maximal stretch with one minimum and at most a cap long.
A record goes to the partition its minimum hashes to.
Pass 2 counts each partition separately and sets edge bits from each k-mer's neighbouring bases.

* `slide_min` proves that the rolling minimum, with its rescans, is the minimum of each window.
  The rescan leaves a stale index when every value in the window is the all-ones word, and the value is still right then.
* By `minOrd_rc`, a k-mer and its reverse complement have the same minimum, so the partition is a function of the canonical k-mer.
* With `records_flatten` and `records_key`, the records hold each k-mer occurrence exactly once, each under its own minimum.
* `count_partition` combines them, so a partition's count of a canonical k-mer is its full count in the input, or 0 when the k-mer belongs to another partition.
* `bits_iff_obs` states that the edge bits are exactly the edges whose (k+1)-mer occurs in the reads, on either strand.
  `obs_symmetric` and `trim_symmetric` give the `hsym` hypothesis of the `Dbg` theorems, before and after trimming to solid neighbours.

The record format, the 2-bit packing and the flank codes are not modelled.
The model assumes pass 2 reads each k-mer occurrence with the bases before and after it in its read, which the flank codes in the record headers supply.
-/

namespace OptikProofs.Count

open OptikProofs.Kmer

/-! ## The rolling window minimum -/

/-- The largest `uint64_t`, the value `best` starts from. -/
def top : ℕ := 2 ^ 64 - 1

/-- `v` is the least of the W values from index `s`. -/
def IsWinMin (h : ℕ → ℕ) (s W v : ℕ) : Prop :=
  (∀ t, s ≤ t → t < s + W → v ≤ h t) ∧ ∃ t, s ≤ t ∧ t < s + W ∧ h t = v

/-- In one comparison of the rescan only a strictly smaller value moves the minimum, so a tie leaves the first in place. -/
def scanStep (h : ℕ → ℕ) (st : ℕ × ℕ) (t : ℕ) : ℕ × ℕ := if h t < st.1 then (h t, t) else st

/-- The rescan over the window from `s`, starting from `best = ~0`. -/
def rescan (h : ℕ → ℕ) (s W bi : ℕ) : ℕ × ℕ := (List.range' s W).foldl (scanStep h) (top, bi)

theorem scan_fold (h : ℕ → ℕ) : ∀ (L : List ℕ) (st : ℕ × ℕ),
    let r := L.foldl (scanStep h) st
    r.1 ≤ st.1 ∧ (∀ t ∈ L, r.1 ≤ h t) ∧ ((r.2 ∈ L ∧ h r.2 = r.1) ∨ r = st)
  | [], st => by simp
  | t :: L, st => by
    obtain ⟨h1, h2, h3⟩ := scan_fold h L (scanStep h st t)
    simp only [List.foldl_cons]
    unfold scanStep at h1 h2 h3 ⊢
    split_ifs at h1 h2 h3 ⊢ with hlt
    · refine ⟨by simp at h1; omega, ?_, ?_⟩
      · intro u hu
        rcases List.mem_cons.mp hu with rfl | hu
        · exact h1
        · exact h2 u hu
      · rcases h3 with ⟨hm, he⟩ | he
        · exact Or.inl ⟨List.mem_cons_of_mem _ hm, he⟩
        · left; rw [he]; simp
    · refine ⟨h1, ?_, ?_⟩
      · intro u hu
        rcases List.mem_cons.mp hu with rfl | hu
        · omega
        · exact h2 u hu
      · rcases h3 with ⟨hm, he⟩ | he
        · exact Or.inl ⟨List.mem_cons_of_mem _ hm, he⟩
        · exact Or.inr he

/-- The invariant of the pair (`best`, `bi`) holds when `best` is the window's minimum and `bi` points at it, unless the minimum is `~0`. -/
def Inv (h : ℕ → ℕ) (s W : ℕ) (st : ℕ × ℕ) : Prop :=
  IsWinMin h s W st.1 ∧ ((s ≤ st.2 ∧ st.2 < s + W ∧ h st.2 = st.1) ∨ st.1 = top)

theorem rescan_inv {h : ℕ → ℕ} (hb : ∀ t, h t ≤ top) {W : ℕ} (hW : 0 < W) (s bi : ℕ) : Inv h s W (rescan h s W bi) := by
  obtain ⟨-, hle, hor⟩ := scan_fold h (List.range' s W) (top, bi)
  simp only [rescan]
  set r := (List.range' s W).foldl (scanStep h) (top, bi)
  have hmem : ∀ t, t ∈ List.range' s W ↔ s ≤ t ∧ t < s + W := fun t => by simp [List.mem_range'_1]
  rcases hor with ⟨hm, he⟩ | he
  · rw [hmem] at hm
    exact ⟨⟨fun t h1 h2 => hle t ((hmem t).mpr ⟨h1, h2⟩), r.2, hm.1, hm.2, he⟩, Or.inl ⟨hm.1, hm.2, he⟩⟩
  · have htop : r.1 = top := by rw [he]
    refine ⟨⟨fun t h1 h2 => hle t ((hmem t).mpr ⟨h1, h2⟩), s, le_rfl, by omega, ?_⟩, Or.inr htop⟩
    have := hle s ((hmem s).mpr ⟨le_rfl, by omega⟩)
    have := hb s
    omega

/-- Advancing one k-mer rescans when the minimum's index has left the window and otherwise compares the new value. -/
def slideStep (h : ℕ → ℕ) (W s : ℕ) (st : ℕ × ℕ) : ℕ × ℕ :=
  if st.2 < s then rescan h s W st.2
  else if h (s + W - 1) < st.1 then (h (s + W - 1), s + W - 1) else st

/-- The state after the k-mer starting at `s` of a run, from any stale index `bi0`. -/
def slide (h : ℕ → ℕ) (W bi0 : ℕ) : ℕ → ℕ × ℕ
  | 0 => rescan h 0 W bi0
  | s + 1 => slideStep h W (s + 1) (slide h W bi0 s)

theorem slideStep_inv {h : ℕ → ℕ} (hb : ∀ t, h t ≤ top) {W : ℕ} (hW : 0 < W) {s : ℕ} {st : ℕ × ℕ}
    (hi : Inv h s W st) : Inv h (s + 1) W (slideStep h W (s + 1) st) := by
  unfold slideStep
  split_ifs with hlt hnew
  · exact rescan_inv hb hW _ _
  · -- The old minimum is still in the window, or it is ~0, so the new minimum is the smaller of it and the new value.
    obtain ⟨⟨hmin, -⟩, hor⟩ := hi
    refine ⟨⟨fun t h1 h2 => ?_, s + W, by omega, by omega, by simp⟩, Or.inl ⟨by omega, by omega, by simp⟩⟩
    by_cases ht : t = s + W
    · subst ht; simp
    · have := hmin t (by omega) (by omega)
      simp only at hnew ⊢
      omega
  · obtain ⟨⟨hmin, -⟩, hor⟩ := hi
    simp only [not_lt] at hlt hnew
    have hW' : s + 1 + W - 1 = s + W := by omega
    rw [hW'] at hnew
    refine ⟨⟨fun t h1 h2 => ?_, ?_⟩, ?_⟩
    · by_cases ht : t = s + W
      · subst ht; exact hnew
      · exact hmin t (by omega) (by omega)
    · rcases hor with ⟨h1, h2, h3⟩ | htop
      · exact ⟨st.2, by omega, by omega, h3⟩
      · refine ⟨s + W, by omega, by omega, ?_⟩
        have := hb (s + W); omega
    · rcases hor with ⟨h1, h2, h3⟩ | htop
      · exact Or.inl ⟨by omega, by omega, h3⟩
      · exact Or.inr htop

/-- The rolling minimum is each window's minimum. -/
theorem slide_min {h : ℕ → ℕ} (hb : ∀ t, h t ≤ top) {W : ℕ} (hW : 0 < W) (bi0 : ℕ) :
    ∀ s, IsWinMin h s W (slide h W bi0 s).1
  | s => (go s).1
where
  go : ∀ s, Inv h s W (slide h W bi0 s)
    | 0 => rescan_inv hb hW 0 bi0
    | s + 1 => slideStep_inv hb hW (go s)

theorem isWinMin_unique {h : ℕ → ℕ} {s W v v' : ℕ} (h1 : IsWinMin h s W v) (h2 : IsWinMin h s W v') : v = v' := by
  obtain ⟨hl1, t1, a1, b1, e1⟩ := h1
  obtain ⟨hl2, t2, a2, b2, e2⟩ := h2
  have := hl1 t2 a2 b2
  have := hl2 t1 a1 b1
  omega

/-! ## The partition is a function of the canonical k-mer -/

section Order

variable (k m : ℕ) (ord : ℕ → ℕ)

/-- The m-mers of a read, indexed by start, as `mer_hashes` holds their order values. -/
theorem mmerAt_window (r : List Base) {s o : ℕ} (ho : o + m ≤ k) :
    Spss.mmerAt m ((r.drop s).take k) o = (r.drop (s + o)).take m := by
  simp only [Spss.mmerAt, List.drop_take, List.take_take, List.drop_drop]
  rw [Nat.add_comm s o, min_eq_left (by omega)]

/-- A k-mer's minimum is its reverse complement's. -/
theorem minOrd_rc {y : List Base} (hy : y.length = k) (hmk : m ≤ k) (v : ℕ) :
    IsWinMin (Spss.ordAt m ord (rc y)) 0 (k - m + 1) v ↔ IsWinMin (Spss.ordAt m ord y) 0 (k - m + 1) v := by
  have hr := fun o (ho : o ≤ k - m) => Spss.ordAt_rc (ord := ord) hy hmk ho
  constructor
  · rintro ⟨hle, t, -, ht, he⟩
    refine ⟨fun o _ ho => ?_, k - m - t, by omega, by omega, ?_⟩
    · have := hle (k - m - o) (by omega) (by omega)
      rw [hr _ (by omega), show k - m - (k - m - o) = o by omega] at this
      exact this
    · rw [← he, hr _ (by omega)]
  · rintro ⟨hle, t, -, ht, he⟩
    refine ⟨fun o _ ho => ?_, k - m - t, by omega, by omega, ?_⟩
    · rw [hr _ (by omega)]; exact hle _ (by omega) (by omega)
    · rw [hr _ (by omega), show k - m - (k - m - t) = t by omega, he]

end Order

/-! ## Records -/

section Records

variable {α : Type*} (key : α → ℕ) (full : List α → Bool)

/-- The emit loop over a run's k-mers extends the open record while the minimum is unchanged and closes it at the cap. -/
def records : List α → List α → List (List α)
  | cur, [] => if cur = [] then [] else [cur]
  | cur, x :: xs =>
    let (out, cur') :=
      match cur.head? with
      | some y => if key y = key x then ([], cur ++ [x]) else ([cur], [x])
      | none => ([], [x])
    if full cur' then out ++ [cur'] ++ records [] xs else out ++ records cur' xs

/-- Every record is non-empty and all its k-mers share one minimum. -/
def KeyOk (g : List α) : Prop := g ≠ [] ∧ ∀ x ∈ g, ∀ y ∈ g, key x = key y

theorem records_spec : ∀ (xs cur : List α), (cur = [] ∨ KeyOk key cur) →
    (records key full cur xs).flatten = cur ++ xs ∧ ∀ g ∈ records key full cur xs, KeyOk key g
  | [], cur, hc => by
    unfold records
    split_ifs with h
    · simp [h]
    · simp only [List.flatten_cons, List.flatten_nil, List.append_nil, List.mem_singleton, forall_eq, true_and]
      rcases hc with h' | h'
      · exact absurd h' h
      · exact h'
  | x :: xs, cur, hc => by
    -- The record being extended or opened satisfies KeyOk, and so does any record it closes.
    have hsingle : KeyOk key [x] := ⟨by simp, by simp⟩
    unfold records
    cases hh : cur.head? with
    | none =>
      have : cur = [] := List.head?_eq_none_iff.mp hh
      subst this
      simp only
      split_ifs
      · obtain ⟨h1, h2⟩ := records_spec xs [] (Or.inl rfl)
        refine ⟨by simp [h1], ?_⟩
        intro g hg
        simp only [List.nil_append, List.singleton_append, List.mem_cons] at hg
        rcases hg with rfl | hg
        · exact hsingle
        · exact h2 g hg
      · obtain ⟨h1, h2⟩ := records_spec xs [x] (Or.inr hsingle)
        exact ⟨by simpa using h1, by simpa using h2⟩
    | some y =>
      have hy : y ∈ cur := List.mem_of_mem_head? hh
      have hcur : KeyOk key cur := by
        rcases hc with h | h
        · subst h; simp at hh
        · exact h
      simp only
      split_ifs with hk hf hf
      · -- Extended and full.
        have hok : KeyOk key (cur ++ [x]) := by
          refine ⟨by simp, fun a ha b hb => ?_⟩
          have hx : ∀ c ∈ cur ++ [x], key c = key y := by
            intro c hc'
            rcases List.mem_append.mp hc' with hc' | hc'
            · exact hcur.2 c hc' y hy
            · simp at hc'; rw [hc', hk]
          rw [hx a ha, hx b hb]
        obtain ⟨h1, h2⟩ := records_spec xs [] (Or.inl rfl)
        refine ⟨by simp [h1], ?_⟩
        intro g hg
        simp only [List.nil_append, List.singleton_append, List.mem_cons] at hg
        rcases hg with rfl | hg
        · exact hok
        · exact h2 g hg
      · have hok : KeyOk key (cur ++ [x]) := by
          refine ⟨by simp, fun a ha b hb => ?_⟩
          have hx : ∀ c ∈ cur ++ [x], key c = key y := by
            intro c hc'
            rcases List.mem_append.mp hc' with hc' | hc'
            · exact hcur.2 c hc' y hy
            · simp at hc'; rw [hc', hk]
          rw [hx a ha, hx b hb]
        obtain ⟨h1, h2⟩ := records_spec xs (cur ++ [x]) (Or.inr hok)
        exact ⟨by simpa using h1, by simpa using h2⟩
      · -- A new minimum closes the open record, and the new one is full at once.
        obtain ⟨h1, h2⟩ := records_spec xs [] (Or.inl rfl)
        refine ⟨by simp [h1], ?_⟩
        intro g hg
        simp only [List.cons_append, List.nil_append, List.mem_cons] at hg
        rcases hg with rfl | rfl | hg
        · exact hcur
        · exact hsingle
        · exact h2 g hg
      · obtain ⟨h1, h2⟩ := records_spec xs [x] (Or.inr hsingle)
        refine ⟨by simp [h1], ?_⟩
        intro g hg
        simp only [List.singleton_append, List.mem_cons] at hg
        rcases hg with rfl | hg
        · exact hcur
        · exact h2 g hg

/-- The records of a run hold its k-mer occurrences once each, in order. -/
theorem records_flatten (xs : List α) : (records key full [] xs).flatten = xs := by
  simpa using (records_spec key full xs [] (Or.inl rfl)).1

/-- Each record's k-mers share its minimum. -/
theorem records_key (xs : List α) : ∀ g ∈ records key full [] xs, KeyOk key g :=
  (records_spec key full xs [] (Or.inl rfl)).2

/-- Sending each record to the partition of its first k-mer's key is sending each k-mer to the partition of its own key. -/
theorem flatten_filter (P : ℕ → Bool) : ∀ gs : List (List α), (∀ g ∈ gs, KeyOk key g) →
    (gs.filter fun g => (g.head?.map (P ∘ key)).getD false).flatten = gs.flatten.filter (P ∘ key)
  | [], _ => rfl
  | g :: gs, h => by
    obtain ⟨hne, hk⟩ := h g (List.mem_cons_self ..)
    have ih := flatten_filter P gs (fun g' hg' => h g' (List.mem_cons_of_mem _ hg'))
    obtain ⟨y, ys, rfl⟩ := List.exists_cons_of_ne_nil hne
    have hall : ∀ x ∈ y :: ys, (P ∘ key) x = (P ∘ key) y := fun x hx => by
      simp only [Function.comp_apply, hk x hx y (List.mem_cons_self ..)]
    rw [List.flatten_cons, List.filter_append, ← ih]
    by_cases hp : (P ∘ key) y
    · rw [List.filter_cons_of_pos (by simpa using hp), List.flatten_cons]
      congr 1
      exact (List.filter_eq_self.mpr fun x hx => by rw [hall x hx]; exact hp).symm
    · rw [List.filter_cons_of_neg (by simpa using hp)]
      have hnil : (y :: ys).filter (P ∘ key) = [] := List.filter_eq_nil_iff.mpr fun x hx => by rw [hall x hx]; exact hp
      rw [hnil, List.nil_append]

end Records

/-! ## Counting per partition -/

section Key

variable (k m : ℕ) (ord : ℕ → ℕ)

/-- A k-mer's minimum order value over its m-mers, as one rescan of its window computes it. -/
def minKey (y : List Base) : ℕ := (rescan (Spss.ordAt m ord y) 0 (k - m + 1) 0).1

variable {k m ord}

theorem minKey_spec (hb : ∀ a, ord a ≤ top) (y : List Base) :
    IsWinMin (Spss.ordAt m ord y) 0 (k - m + 1) (minKey k m ord y) :=
  (rescan_inv (fun _ => hb _) (by omega) 0 0).1

theorem minKey_rc (hb : ∀ a, ord a ≤ top) (hmk : m ≤ k) {y : List Base} (hy : y.length = k) :
    minKey k m ord (rc y) = minKey k m ord y :=
  isWinMin_unique ((minOrd_rc k m ord hy hmk _).mp (minKey_spec hb (rc y))) (minKey_spec hb y)

/-- k-mers with one canonical form have one minimum, so they go to one partition. -/
theorem minKey_canon (hb : ∀ a, ord a ≤ top) (hmk : m ≤ k) {y z : List Base} (hy : y.length = k) (hz : z.length = k)
    (h : canon y = canon z) : minKey k m ord y = minKey k m ord z := by
  rcases (canon_eq_iff (hy.trans hz.symm)).mp h with rfl | rfl
  · rfl
  · exact (minKey_rc hb hmk hy).symm

/-- The rolling minimum over a run's m-mer values, at the k-mer starting at `s`, is that k-mer's `minKey`. -/
theorem slide_eq_minKey (hb : ∀ a, ord a ≤ top) (hmk : m ≤ k) (r : List Base) (bi0 s : ℕ) :
    (slide (fun t => ord (canon ((r.drop t).take m))) (k - m + 1) bi0 s).1 = minKey k m ord ((r.drop s).take k) := by
  have hshift : ∀ o, o < k - m + 1 → Spss.ordAt m ord ((r.drop s).take k) o = ord (canon ((r.drop (s + o)).take m)) := by
    intro o ho
    simp only [Spss.ordAt]
    rw [mmerAt_window k m r (by omega)]
  obtain ⟨hle, t, hst, ht, he⟩ := slide_min (h := fun t => ord (canon ((r.drop t).take m))) (fun _ => hb _)
    (W := k - m + 1) (by omega) bi0 s
  refine isWinMin_unique ⟨fun o _ ho => ?_, t - s, by omega, by omega, ?_⟩ (minKey_spec hb _)
  · rw [hshift o (by omega)]; exact hle (s + o) (by omega) (by omega)
  · rw [hshift _ (by omega), show s + (t - s) = t by omega]; exact he

end Key

theorem flatten_flatMap {α : Type*} (f : List α → List (List α)) (hf : ∀ r, (f r).flatten = r) :
    ∀ runs : List (List α), (runs.flatMap f).flatten = runs.flatten
  | [] => rfl
  | r :: runs => by rw [List.flatMap_cons, List.flatten_append, hf, flatten_flatMap f hf runs, List.flatten_cons]

/-- Partitioned counting is exact.
Partition p counts the records whose minimum hashes to p, which hold every occurrence of a canonical k-mer that belongs to p and none of any other. -/
theorem count_partition (keyOf : List Base → ℕ) (full : List (List Base) → Bool) (part : ℕ → ℕ)
    (runs : List (List (List Base))) (c : List Base) (p : ℕ)
    (hkey : ∀ y ∈ runs.flatten, canon y = canon c → keyOf y = keyOf c) :
    let P : ℕ → Bool := fun v => decide (part v = p)
    (((runs.flatMap (records keyOf full [])).filter fun g => (g.head?.map (P ∘ keyOf)).getD false).flatten).countP
        (fun y => canon y = canon c) =
      if part (keyOf c) = p then runs.flatten.countP (fun y => canon y = canon c) else 0 := by
  intro P
  have hok : ∀ g ∈ runs.flatMap (records keyOf full []), KeyOk keyOf g := by
    intro g hg
    obtain ⟨r, -, hr⟩ := List.mem_flatMap.mp hg
    exact records_key keyOf full r g hr
  rw [flatten_filter keyOf P _ hok, flatten_flatMap _ (records_flatten keyOf full), List.countP_filter]
  split_ifs with hp
  · apply List.countP_congr
    intro y hy
    simp only [Bool.and_eq_true, decide_eq_true_eq, Function.comp_apply, P]
    constructor
    · exact fun h => h.1
    · exact fun h => ⟨h, by rw [hkey y hy h]; exact hp⟩
  · rw [List.countP_eq_zero]
    intro y hy
    simp only [Bool.and_eq_true, decide_eq_true_eq, Function.comp_apply, P, not_and]
    intro h
    rw [hkey y hy h]; exact hp

/-! ## Edge bits -/

/-- The (k+1)-mer a bit stands for at node c, which is c followed by b for out-bit b and a followed by c for in-bit a. -/
def word (c : List Base) : Bool × Base → List Base
  | (true, b) => c ++ [b]
  | (false, a) => a :: c

/-- The node a k-mer's bits go to, its canonical orientation. -/
def node (y : List Base) : List Base := if pack y ≤ pack (rc y) then y else rc y

/-- The bits one occurrence sets, each an out flag and a base, as pass 2 builds `in_b` and `out_b`.
On the reverse strand the neighbours swap sides and are complemented, and a palindrome takes both branches. -/
def occBits (y : List Base) (prev next : Option Base) : List (Bool × Base) :=
  (if pack y ≤ pack (rc y) then next.toList.map (true, ·) ++ prev.toList.map (false, ·) else []) ++
  (if pack (rc y) ≤ pack y then next.toList.map (fun b => (false, comp b)) ++ prev.toList.map (fun a => (true, comp a))
   else [])

variable (k : ℕ)

/-- The occurrence of the k-mer at `s` in a read, with the bases before and after it if the read has them. -/
def occAt (r : List Base) (s : ℕ) : List Base × Option Base × Option Base :=
  ((r.drop s).take k, (if s = 0 then none else r[s - 1]?), r[s + k]?)

/-- Bit `bt` is set at node c when some occurrence whose node is c sets it. -/
def Bit (reads : List (List Base)) (c : List Base) (bt : Bool × Base) : Prop :=
  ∃ r ∈ reads, ∃ s, s + k ≤ r.length ∧
    node (occAt k r s).1 = c ∧ bt ∈ occBits (occAt k r s).1 (occAt k r s).2.1 (occAt k r s).2.2

/-- The (k+1)-mers of the reads. -/
def Obs (reads : List (List Base)) (w : List Base) : Prop := ∃ r ∈ reads, ∃ s, s + (k + 1) ≤ r.length ∧ (r.drop s).take (k + 1) = w

variable {k}

theorem node_of_le {y : List Base} (h : pack y ≤ pack (rc y)) : node y = y := by simp [node, h]

theorem node_of_ge {y : List Base} (h : pack (rc y) ≤ pack y) : node y = rc y := by
  unfold node
  split_ifs with h'
  · exact pack_injective (x := y) (y := rc y) (by simp) (le_antisymm h' h)
  · rfl

theorem window_succ {r : List Base} {s : ℕ} (h : s + (k + 1) ≤ r.length) :
    (r.drop s).take (k + 1) = (r.drop s).take k ++ [r[s + k]'(by omega)] := by
  rw [List.take_add_one, List.getElem?_drop]
  simp [List.getElem?_eq_getElem (show s + k < r.length by omega)]

theorem window_cons {r : List Base} {s : ℕ} (h : s + (k + 1) ≤ r.length) :
    (r.drop s).take (k + 1) = r[s]'(by omega) :: (r.drop (s + 1)).take k := by
  rw [List.drop_eq_getElem_cons (by omega : s < r.length), List.take_succ_cons]

theorem rc_append_single (x : List Base) (b : Base) : rc (x ++ [b]) = comp b :: rc x := by simp [rc]
theorem rc_cons (a : Base) (x : List Base) : rc (a :: x) = rc x ++ [comp a] := by simp [rc]

theorem mem_occBits {y : List Base} {prev next : Option Base} {bt : Bool × Base} :
    bt ∈ occBits y prev next ↔
      (pack y ≤ pack (rc y) ∧ ((∃ b, next = some b ∧ bt = (true, b)) ∨ (∃ a, prev = some a ∧ bt = (false, a)))) ∨
      (pack (rc y) ≤ pack y ∧ ((∃ b, next = some b ∧ bt = (false, comp b)) ∨ (∃ a, prev = some a ∧ bt = (true, comp a)))) := by
  unfold occBits
  by_cases h1 : pack y ≤ pack (rc y) <;> by_cases h2 : pack (rc y) ≤ pack y <;>
    cases next <;> cases prev <;> simp [h1, h2, or_assoc]

theorem getElem_of_some {r : List Base} {i : ℕ} {b : Base} (h : r[i]? = some b) : ∃ hi : i < r.length, r[i] = b := by
  rw [List.getElem?_eq_some_iff] at h; exact h

/-- An edge bit is set exactly when its (k+1)-mer occurs in the reads, on one strand or the other. -/
theorem bits_iff_obs (reads : List (List Base)) {c : List Base} (hc : c.length = k) (hn : node c = c)
    (bt : Bool × Base) : Bit k reads c bt ↔ Obs k reads (word c bt) ∨ Obs k reads (rc (word c bt)) := by
  have hle : pack c ≤ pack (rc c) := by
    unfold node at hn
    split_ifs at hn with h
    · exact h
    · rw [hn]
  constructor
  · rintro ⟨r, hr, s, hs, hnode, hbt⟩
    simp only [occAt] at hnode hbt
    rcases mem_occBits.mp hbt with ⟨h1, ⟨b, hb, rfl⟩ | ⟨a, ha, rfl⟩⟩ | ⟨h2, ⟨b, hb, rfl⟩ | ⟨a, ha, rfl⟩⟩
    · -- On the forward strand the next base comes from the window at s.
      rw [node_of_le h1] at hnode; subst hnode
      obtain ⟨hi, rfl⟩ := getElem_of_some hb
      exact Or.inl ⟨r, hr, s, by omega, by rw [window_succ (by omega)]; rfl⟩
    · -- The previous base on the forward strand comes from the window at s - 1.
      rw [node_of_le h1] at hnode; subst hnode
      split_ifs at ha with h0
      obtain ⟨hi, rfl⟩ := getElem_of_some ha
      refine Or.inl ⟨r, hr, s - 1, by omega, ?_⟩
      rw [window_cons (by omega), show s - 1 + 1 = s by omega]; rfl
    · rw [node_of_ge h2] at hnode; subst hnode
      obtain ⟨hi, rfl⟩ := getElem_of_some hb
      refine Or.inr ⟨r, hr, s, by omega, ?_⟩
      rw [window_succ (by omega)]
      simp [word, rc_cons]
    · rw [node_of_ge h2] at hnode; subst hnode
      split_ifs at ha with h0
      obtain ⟨hi, rfl⟩ := getElem_of_some ha
      refine Or.inr ⟨r, hr, s - 1, by omega, ?_⟩
      rw [window_cons (by omega), show s - 1 + 1 = s by omega]
      simp [word, rc_append_single]
  · intro h
    obtain ⟨o, b⟩ := bt
    have hlen : ∀ (r : List Base) s, s + (k + 1) ≤ r.length → ((r.drop s).take k).length = k := by
      intro r s hs; simp; omega
    rcases h with ⟨r, hr, s, hs, hw⟩ | ⟨r, hr, s, hs, hw⟩ <;> cases o
    · -- For b followed by c, the occurrence at s + 1 has b before it.
      rw [window_cons hs] at hw
      simp only [word, List.cons.injEq] at hw
      obtain ⟨rfl, hw⟩ := hw
      refine ⟨r, hr, s + 1, by omega, ?_, ?_⟩
      · simp only [occAt]; rw [hw]; exact hn
      · simp only [occAt]; rw [hw]
        exact mem_occBits.mpr (Or.inl ⟨hle, Or.inr ⟨_, by simp [List.getElem?_eq_getElem (show s < r.length by omega)], rfl⟩⟩)
    · -- When c is followed by b, the occurrence at s has b after it.
      rw [window_succ hs] at hw
      simp only [word] at hw
      obtain ⟨hw, hb⟩ := List.append_inj' hw rfl
      simp only [List.cons.injEq, and_true] at hb
      refine ⟨r, hr, s, by omega, ?_, ?_⟩
      · simp only [occAt]; rw [hw]; exact hn
      · simp only [occAt]; rw [hw]
        exact mem_occBits.mpr (Or.inl ⟨hle, Or.inl ⟨_, by rw [List.getElem?_eq_getElem (show s + k < r.length by omega), hb], rfl⟩⟩)
    · -- The other strand of b followed by c is rc c followed by comp b, which the occurrence at s reads.
      rw [window_succ hs] at hw
      simp only [word, rc_cons] at hw
      obtain ⟨hw, hb⟩ := List.append_inj' hw rfl
      simp only [List.cons.injEq, and_true] at hb
      have hge : pack (rc (rc c)) ≤ pack (rc c) := by rw [rc_rc]; exact hle
      refine ⟨r, hr, s, by omega, ?_, ?_⟩
      · simp only [occAt]; rw [hw, node_of_ge hge, rc_rc]
      · simp only [occAt]; rw [hw]
        refine mem_occBits.mpr (Or.inr ⟨hge, Or.inl ⟨comp b, ?_, by simp⟩⟩)
        rw [List.getElem?_eq_getElem (show s + k < r.length by omega), hb]
    · -- The other strand of c followed by b is comp b followed by rc c, which the occurrence at s + 1 reads.
      rw [window_cons hs] at hw
      simp only [word, rc_append_single, List.cons.injEq] at hw
      obtain ⟨hb, hw⟩ := hw
      have hge : pack (rc (rc c)) ≤ pack (rc c) := by rw [rc_rc]; exact hle
      refine ⟨r, hr, s + 1, by omega, ?_, ?_⟩
      · simp only [occAt]; rw [hw, node_of_ge hge, rc_rc]
      · simp only [occAt]; rw [hw]
        refine mem_occBits.mpr (Or.inr ⟨hge, Or.inr ⟨comp b, ?_, by simp⟩⟩)
        simp [List.getElem?_eq_getElem (show s < r.length by omega), hb]

/-! ## The edge relation -/

variable (k)

/-- The last base of a k-mer. -/
def lastB (X : List Base) : Base := X.getLastD 0

/-- The edge from x to X as the graph reads it at x's node.
It is an out-bit when x is in canonical orientation, and otherwise the complemented in-bit of the other strand. -/
def CodeEdge (reads : List (List Base)) (x X : List Base) : Prop :=
  Dbg.Overlap k x X ∧
    Bit k reads (node x) (if pack x ≤ pack (rc x) then (true, lastB X) else (false, comp (lastB X)))

/-- The edge from x to X when their (k+1)-mer occurs in the reads on either strand, the relation the Dbg vectors use. -/
def ObsEdge (reads : List (List Base)) (x X : List Base) : Prop :=
  Dbg.Overlap k x X ∧ (Obs k reads (Dbg.rightWord k x X) ∨ Obs k reads (rc (Dbg.rightWord k x X)))

variable {k}

theorem rightWord_last (hk : 1 ≤ k) {x X : List Base} (hX : X.length = k) : Dbg.rightWord k x X = x ++ [lastB X] := by
  unfold Dbg.rightWord lastB
  congr 1
  obtain ⟨Y, b, rfl⟩ : ∃ Y b, X = Y ++ [b] := ⟨X.dropLast, X.getLast (by intro h; simp [h] at hX; omega),
    (List.dropLast_append_getLast _).symm⟩
  simp at hX
  rw [List.getLastD_eq_getLast?, List.getLast?_append]
  simp [show k - 1 = Y.length by omega]

/-- The bits are the observed edges. -/
theorem codeEdge_iff (hk : 1 ≤ k) (reads : List (List Base)) (x X : List Base) :
    CodeEdge k reads x X ↔ ObsEdge k reads x X := by
  unfold CodeEdge ObsEdge
  refine and_congr_right fun ho => ?_
  rw [rightWord_last hk ho.2.1]
  split_ifs with h
  · rw [node_of_le h, bits_iff_obs reads ho.1 (node_of_le h)]
    rfl
  · have h' : pack (rc x) ≤ pack x := le_of_lt (not_le.mp h)
    have hn : node (rc x) = rc x := node_of_le (by rw [rc_rc]; exact h')
    rw [node_of_ge h', bits_iff_obs reads (by simp [ho.1]) hn]
    simp only [word, ← rc_append_single, rc_rc]
    exact or_comm

theorem overlap_rc {x X : List Base} (h : Dbg.Overlap k x X) : Dbg.Overlap k (rc X) (rc x) := by
  obtain ⟨hx, hX, ht⟩ := h
  refine ⟨by simp [hX], by simp [hx], ?_⟩
  simp only [rc, List.tail_reverse, List.dropLast_reverse, ← List.map_dropLast]
  rw [← ht]
  simp [List.map_tail]

theorem rightWord_rc (hk : 1 ≤ k) {x X : List Base} (h : Dbg.Overlap k x X) :
    Dbg.rightWord k (rc X) (rc x) = rc (Dbg.rightWord k x X) := by
  rw [Dbg.rightWord_eq_leftWord h, Dbg.rightWord, Dbg.drop_rc, h.1, show k - (k - 1) = 1 by omega]
  simp [rc, Dbg.leftWord, List.map_take]

/-- The observed edges are strand-symmetric, the `hsym` hypothesis of `Dbg.compacts_rc` and `Dbg.unitig_rc`. -/
theorem obs_symmetric (hk : 1 ≤ k) (reads : List (List Base)) (x X : List Base) :
    ObsEdge k reads x X ↔ ObsEdge k reads (rc X) (rc x) := by
  have one : ∀ x X, ObsEdge k reads x X → ObsEdge k reads (rc X) (rc x) := by
    rintro x X ⟨ho, hw⟩
    refine ⟨overlap_rc ho, ?_⟩
    rw [rightWord_rc hk ho, rc_rc]
    exact hw.symm
  exact ⟨one x X, fun h => by simpa using one _ _ h⟩

theorem code_symmetric (hk : 1 ≤ k) (reads : List (List Base)) (x X : List Base) :
    CodeEdge k reads x X ↔ CodeEdge k reads (rc X) (rc x) := by
  rw [codeEdge_iff hk, codeEdge_iff hk, obs_symmetric hk]

variable (k) in
/-- With min_count above 1, the pass-2 trim and the solid-node filter remove every edge with a non-solid end. -/
def TrimEdge (solid : ℕ → Prop) (reads : List (List Base)) (x X : List Base) : Prop :=
  CodeEdge k reads x X ∧ solid (canon x) ∧ solid (canon X)

theorem trim_symmetric (hk : 1 ≤ k) (solid : ℕ → Prop) (reads : List (List Base)) (x X : List Base) :
    TrimEdge k solid reads x X ↔ TrimEdge k solid reads (rc X) (rc x) := by
  unfold TrimEdge
  rw [code_symmetric hk, canon_rc, canon_rc]
  tauto

end OptikProofs.Count
