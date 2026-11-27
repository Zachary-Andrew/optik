import OptikProofs.Kmer

/-!
# No FM match crosses a unitig separator

`FmIndex` in `include/fm_index.hpp` indexes every unitig followed by the separator X.

A symbol is a base or X, written `none`.
A query is ACGT only, so it contains no X.
`occs_sep` shows that its occurrences in a text with an X in it are the ones on either side, so no match crosses the X.
By `count_text`, a pattern occurs in the indexed text as many times as it occurs in all the unitigs together.
Then `occurs_text` states that a k-mer is found exactly when some unitig contains it.
With the unitigs a partition of the k-mer set (`Dbg.unitig_starts`, and the vectors' unitig checks), the index answers exactly the k-mer set.
-/

namespace OptikProofs.Fm

open OptikProofs.Kmer

/-- A base, or the separator X as `none`. -/
abbrev Sym := Option Base

/-- The start positions where `y` occurs in `t`. -/
def occs (t y : List Sym) : List ℕ := (List.range t.length).filter fun p => decide (y <+: t.drop p)

/-- The indexed text, each unitig followed by X. -/
def text (us : List (List Base)) : List Sym := us.flatMap fun u => u.map some ++ [none]

/-- A pattern without X is a prefix of `c ++ X :: b` only if it is a prefix of `c`. -/
theorem prefix_sep {y : List Sym} (hy : none ∉ y) : ∀ (c b : List Sym), y <+: c ++ none :: b ↔ y <+: c
  | [], b => by
    constructor
    · intro h
      cases y with
      | nil => exact List.nil_prefix
      | cons a y =>
        obtain ⟨t, ht⟩ := h
        simp only [List.nil_append, List.cons_append, List.cons.injEq] at ht
        exact absurd (ht.1 ▸ List.mem_cons_self ..) hy
    · intro h
      rw [List.prefix_nil] at h
      subst h
      exact List.nil_prefix
  | a :: c, b => by
    cases y with
    | nil => simp
    | cons y0 y =>
      have hy' : none ∉ y := fun h => hy (List.mem_cons_of_mem _ h)
      simp only [List.cons_append, List.cons_prefix_cons, prefix_sep hy' c b]

theorem mem_occs {t y : List Sym} {p : ℕ} : p ∈ occs t y ↔ p < t.length ∧ y <+: t.drop p := by
  simp [occs]

theorem countP_range_add (f : ℕ → Bool) (m k : ℕ) :
    (List.range (m + k)).countP f = (List.range m).countP f + (List.range k).countP fun q => f (m + q) := by
  rw [List.range_add, List.countP_append, List.countP_map]
  rfl

theorem length_occs (t y : List Sym) : (occs t y).length = (List.range t.length).countP fun p => decide (y <+: t.drop p) := by
  rw [occs, List.countP_eq_length_filter]

/-- An X splits the occurrences of a pattern without X into those in `a` and those in `b`, with none across it. -/
theorem occs_sep {y : List Sym} (hy : none ∉ y) (hne : y ≠ []) (a b : List Sym) :
    (occs (a ++ none :: b) y).length = (occs a y).length + (occs b y).length := by
  rw [length_occs, length_occs, length_occs]
  have hlen : (a ++ none :: b).length = a.length + (1 + b.length) := by simp; omega
  rw [hlen, countP_range_add, countP_range_add]
  have h1 : ((List.range a.length).countP fun p => decide (y <+: (a ++ none :: b).drop p)) =
      (List.range a.length).countP fun p => decide (y <+: a.drop p) := by
    apply List.countP_congr
    intro p hp
    rw [List.mem_range] at hp
    rw [List.drop_append_of_le_length (by omega)]
    simp only [decide_eq_true_eq]
    exact prefix_sep hy _ _
  have h2 : ((List.range 1).countP fun q => decide (y <+: (a ++ none :: b).drop (a.length + q))) = 0 := by
    have : ¬ y <+: none :: b := by
      intro h
      have := (prefix_sep hy [] b).mp h
      rw [List.prefix_nil] at this
      exact hne this
    simp [this]
  have h3 : ((List.range b.length).countP fun q => decide (y <+: (a ++ none :: b).drop (a.length + (1 + q)))) =
      (List.range b.length).countP fun p => decide (y <+: b.drop p) := by
    apply List.countP_congr
    intro p _
    have hd : (a ++ none :: b).drop (a.length + 1) = b := by rw [← List.drop_drop, List.drop_left' rfl]; rfl
    rw [show a.length + (1 + p) = a.length + 1 + p by omega, ← List.drop_drop, hd]
  rw [h1, h2, h3, Nat.zero_add]

theorem text_cons (u : List Base) (us : List (List Base)) : text (u :: us) = u.map some ++ none :: text us := by
  simp [text]

/-- A pattern without X occurs in the indexed text as often as in all the unitigs together. -/
theorem count_text {y : List Sym} (hy : none ∉ y) (hne : y ≠ []) :
    ∀ us : List (List Base), (occs (text us) y).length = (us.map fun u => (occs (u.map some) y).length).sum
  | [] => by simp [text, occs]
  | u :: us => by rw [text_cons, occs_sep hy hne, count_text hy hne us, List.map_cons, List.sum_cons]

theorem prefix_map_some : ∀ (x l : List Base), x.map some <+: l.map some ↔ x <+: l
  | [], _ => by simp
  | _ :: _, [] => by simp
  | a :: x, b :: l => by simp [List.cons_prefix_cons, prefix_map_some x l]

/-- A k-mer is found in the indexed text exactly when some unitig contains it. -/
theorem occurs_text {x : List Base} (hne : x ≠ []) (us : List (List Base)) :
    (∃ p, p < (text us).length ∧ x.map some <+: (text us).drop p) ↔ ∃ u ∈ us, ∃ p, p < u.length ∧ x <+: u.drop p := by
  have hy : none ∉ x.map some := by simp
  have hne' : x.map some ≠ [] := by simpa using hne
  have hc := count_text hy hne' us
  have hpos : ∀ t, (∃ p, p < t.length ∧ x.map some <+: t.drop p) ↔ 0 < (occs t (x.map some)).length := by
    intro t
    rw [List.length_pos_iff_exists_mem]
    simp only [mem_occs]
  rw [hpos, hc]
  have hu : ∀ u : List Base, 0 < (occs (u.map some) (x.map some)).length ↔ ∃ p, p < u.length ∧ x <+: u.drop p := by
    intro u
    rw [← hpos]
    simp [← List.map_drop, prefix_map_some]
  have hsum : ∀ vs : List (List Base), 0 < (vs.map fun u => (occs (u.map some) (x.map some)).length).sum ↔
      ∃ u ∈ vs, 0 < (occs (u.map some) (x.map some)).length := by
    intro vs
    induction vs with
    | nil => simp
    | cons u vs ih =>
      rw [List.map_cons, List.sum_cons]
      constructor
      · intro h
        by_cases h0 : 0 < (occs (u.map some) (x.map some)).length
        · exact ⟨u, List.mem_cons_self .., h0⟩
        · obtain ⟨w, hw, hw'⟩ := ih.mp (by omega)
          exact ⟨w, List.mem_cons_of_mem _ hw, hw'⟩
      · rintro ⟨w, hw, hw'⟩
        rcases List.mem_cons.mp hw with rfl | hw
        · omega
        · have := ih.mpr ⟨w, hw, hw'⟩
          omega
  rw [hsum]
  exact exists_congr fun u => and_congr_right fun _ => hu u

end OptikProofs.Fm
