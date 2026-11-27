import Mathlib.Data.Finset.Card
import Mathlib.Order.Interval.Finset.Nat

/-!
# Elias-Fano keys and ranks are exact

Models `DbgGraphEF::key` and `DbgGraphEF::rank` in `include/dbg_graph_ef.hpp`, which store the graph's sorted canonical k-mers, and `EfSeq::operator[]` in `include/phobic_mphf.hpp`.

Key i is split into its high part `v i >> l` and its low l bits.
The high parts form a unary bitvector with a one at position `hi i + i` for every key.
`key i` reads the i-th one with select1 and takes away i, and `rank c` finds the bucket of c's high part with select0 and scans it for c's low bits.

* `key_eq` shows that `key i` is the i-th key, for any keys in non-decreasing order, which covers `EfSeq` too.
* `rank_some` and `rank_complete` prove that `rank c` returns i exactly when key i is c, for keys in strictly increasing order.
  `dbg_graph_ef.hpp` needs this for "membership is exact" and "node i is the i-th smallest k-mer".

select1 and select0 are taken by their definition, the position of the i-th one or zero.
The code finds them from a sample every 256 ones or zeros and a popcount walk over words, which is not modelled here.
-/

namespace OptikProofs.Ef

open Finset

variable (n l : ℕ) (v : ℕ → ℕ)

/-- The high part of key j. -/
def hi (j : ℕ) : ℕ := v j / 2 ^ l

/-- The low part of key j, `low(j)` in the code. -/
def lo (j : ℕ) : ℕ := v j % 2 ^ l

/-- Bit p of the high bitvector is set. -/
def H (p : ℕ) : Prop := ∃ j < n, hi l v j + j = p

instance (p : ℕ) : Decidable (H n l v p) := by unfold H; infer_instance

/-- The ones below position p. -/
def ones (p : ℕ) : ℕ := ((range n).filter fun j => hi l v j + j < p).card

/-- p is the select1 answer for i, the position of the one with i ones before it. -/
def Sel1 (i p : ℕ) : Prop := H n l v p ∧ ones n l v p = i

/-- The zero at p, the select0 answer for t, has t zeros before it. -/
def Sel0 (t p : ℕ) : Prop := ¬H n l v p ∧ p - ones n l v p = t

variable {n l v}

theorem hi_mono (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {i j : ℕ} (hij : i ≤ j) (hj : j < n) : hi l v i ≤ hi l v j :=
  Nat.div_le_div_right (hs i j hij hj)

/-- Ones sit at strictly increasing positions. -/
theorem pos_lt_iff (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {i j : ℕ} (hi' : i < n) (hj : j < n) :
    hi l v j + j < hi l v i + i ↔ j < i := by
  constructor
  · intro h
    by_contra hc
    have := hi_mono (l := l) hs (show i ≤ j by omega) hj
    omega
  · intro h
    have := hi_mono (l := l) hs (show j ≤ i by omega) hi'
    omega

theorem ones_pos (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {i : ℕ} (hi' : i < n) : ones n l v (hi l v i + i) = i := by
  unfold ones
  rw [show ((range n).filter fun j => hi l v j + j < hi l v i + i) = range i by
    ext j
    simp only [mem_filter, mem_range]
    constructor
    · rintro ⟨hj, h⟩; exact (pos_lt_iff (l := l) hs hi' hj).mp h
    · intro h; exact ⟨by omega, (pos_lt_iff hs hi' (by omega)).mpr h⟩]
  simp

/-- `key(i)` takes the i-th one's position less i as the high part and ORs the low part below it. -/
theorem key_eq (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {i p : ℕ} (hi' : i < n) (h : Sel1 n l v i p) :
    (p - i) <<< l ||| lo l v i = v i := by
  obtain ⟨⟨j, hj, rfl⟩, hcnt⟩ := h
  rw [ones_pos hs hj] at hcnt
  subst hcnt
  rw [Nat.add_sub_cancel, show lo l v j = v j % 2 ^ l from rfl,
    ← Nat.shiftLeft_add_eq_or_of_lt (Nat.mod_lt _ (Nat.two_pow_pos l)), Nat.shiftLeft_eq]
  unfold hi
  exact Nat.div_add_mod' (v j) (2 ^ l)

/-! ## rank -/

/-- The keys whose high part is below h, which with sorted keys are the first ones. -/
def below (h : ℕ) : ℕ := ((range n).filter fun j => hi l v j < h).card

theorem below_iff (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) (h : ℕ) {j : ℕ} (hj : j < n) :
    hi l v j < h ↔ j < below (n := n) (l := l) (v := v) h := by
  -- The keys below h form a prefix, because high parts never decrease.
  have hpre : ((range n).filter fun j => hi l v j < h) = range ((range n).filter fun j => hi l v j < h).card := by
    set S := (range n).filter fun j => hi l v j < h
    have hdown : ∀ a b, a ≤ b → b ∈ S → a ∈ S := by
      intro a b hab hb
      simp only [S, mem_filter, mem_range] at hb ⊢
      exact ⟨by omega, lt_of_le_of_lt (hi_mono (l := l) hs hab hb.1) hb.2⟩
    -- A downward-closed set of naturals is the range of its size.
    ext a
    simp only [mem_range]
    constructor
    · intro ha
      by_contra hc
      have : range (a + 1) ⊆ S := fun b hb => hdown b a (by simp at hb; omega) ha
      have := card_le_card this
      simp at this
      omega
    · intro ha
      by_contra hn
      have : S ⊆ range a := by
        intro b hb
        simp only [mem_range]
        by_contra hc
        exact hn (hdown a b (by omega) hb)
      have := card_le_card this
      simp at this
      omega
  have hm : j ∈ (range n).filter (fun j => hi l v j < h) ↔ hi l v j < h := by simp [hj]
  rw [← hm, hpre, mem_range]
  rfl

theorem below_le (h : ℕ) : below (n := n) (l := l) (v := v) h ≤ n := by
  unfold below
  exact (card_le_card (filter_subset _ _)).trans (by simp)

/-- The zero before bucket h ≥ 1 sits at h - 1 plus the keys below h, so the bucket starts at index `below h`. -/
theorem sel0_bucket (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {h p : ℕ} (h1 : 1 ≤ h) (hp : Sel0 n l v (h - 1) p) :
    p + 1 - h = below (n := n) (l := l) (v := v) h := by
  set b := below (n := n) (l := l) (v := v) h
  have hb := below_le (n := n) (l := l) (v := v) h
  set q := h - 1 + b
  -- q is a zero with h - 1 zeros before it.
  have hq_ones : ones n l v q = b := by
    unfold ones
    rw [show ((range n).filter fun j => hi l v j + j < q) = range b by
      ext j
      simp only [mem_filter, mem_range]
      constructor
      · rintro ⟨hj, hlt⟩
        by_contra hc
        have := (below_iff hs h hj).not.mpr (by omega)
        omega
      · intro hj
        have := (below_iff hs h (show j < n by omega)).mpr hj
        exact ⟨by omega, by omega⟩]
    simp
  have hq_zero : ¬H n l v q := by
    rintro ⟨j, hj, he⟩
    by_cases hc : j < b
    · have := (below_iff hs h hj).mpr hc; omega
    · have := (below_iff hs h hj).not.mpr hc; omega
  -- Two zeros with the same number of zeros before them are one position.
  have hones_succ : ∀ x, ones n l v (x + 1) = ones n l v x + if H n l v x then 1 else 0 := by
    intro x
    unfold ones
    rw [show ((range n).filter fun j => hi l v j + j < x + 1) =
        ((range n).filter fun j => hi l v j + j < x) ∪ ((range n).filter fun j => hi l v j + j = x) by
      ext j; simp only [mem_filter, mem_union, mem_range]; omega]
    rw [card_union_of_disjoint (disjoint_filter.mpr fun j _ h1 h2 => by omega)]
    congr 1
    by_cases hx : H n l v x
    · rw [if_pos hx]
      obtain ⟨i, hi', rfl⟩ := hx
      rw [card_eq_one]
      refine ⟨i, ?_⟩
      ext j
      simp only [mem_filter, mem_range, mem_singleton]
      constructor
      · rintro ⟨hj, he⟩
        rcases lt_trichotomy j i with hji | rfl | hji
        · have := (pos_lt_iff (l := l) hs hi' hj).mpr hji; omega
        · rfl
        · have := (pos_lt_iff (l := l) hs hj hi').mpr hji; omega
      · rintro rfl; exact ⟨hi', rfl⟩
    · rw [if_neg hx, card_eq_zero, filter_eq_empty_iff]
      intro j hj he
      exact hx ⟨j, by simpa using hj, he⟩
  have hzero_mono : ∀ x y, x ≤ y → x - ones n l v x ≤ y - ones n l v y := by
    intro x y hxy
    induction y, hxy using Nat.le_induction with
    | base => exact le_rfl
    | succ y _ ih =>
      have := hones_succ y
      have hle : ones n l v y ≤ y := by
        unfold ones
        calc _ ≤ (range y).card := card_le_card (fun j hj => by
                simp only [mem_filter, mem_range] at hj ⊢; omega)
          _ = y := card_range y
      by_cases hy : H n l v y
      · rw [if_pos hy] at this; omega
      · rw [if_neg hy] at this; omega
  have hzero_step : ∀ x, ¬H n l v x → (x + 1) - ones n l v (x + 1) = x - ones n l v x + 1 := by
    intro x hx
    have := hones_succ x
    have hle : ones n l v x ≤ x := by
      unfold ones
      calc _ ≤ (range x).card := card_le_card (fun j hj => by
              simp only [mem_filter, mem_range] at hj ⊢; omega)
        _ = x := card_range x
    simp only [hx, if_false] at this
    omega
  obtain ⟨hpz, hpc⟩ := hp
  have hpq : p = q := by
    rcases lt_trichotomy p q with hlt | heq | hgt
    · have h1' := hzero_mono (p + 1) q (by omega)
      rw [hzero_step p hpz] at h1'
      rw [hq_ones] at h1'
      omega
    · exact heq
    · have h1' := hzero_mono (q + 1) p (by omega)
      rw [hzero_step q hq_zero, hq_ones] at h1'
      omega
  omega

/-- The bucket scan, from key index `idx` at position `h + idx`, with `fuel` steps left.
It stops at the first set bit that is not followed, at a matching low part, or at a larger one. -/
def scan (h low : ℕ) : ℕ → ℕ → Option ℕ
  | 0, _ => none
  | fuel + 1, idx =>
    if idx < n ∧ H n l v (h + idx) then
      if lo l v idx = low then some idx else if lo l v idx > low then none else scan h low fuel (idx + 1)
    else none

/-- `rank(c)` given the select0 answer `p0` for bucket h = c >> l, or the bucket start 0 when h = 0. -/
def rank (c : ℕ) (p0 : ℕ) : Option ℕ :=
  let h := c / 2 ^ l
  let pos := if h = 0 then 0 else p0 + 1
  scan (n := n) (l := l) (v := v) h (c % 2 ^ l) (n + 1) (pos - h)

theorem split_key (x : ℕ) : x = x / 2 ^ l * 2 ^ l + x % 2 ^ l := by
  have := Nat.div_add_mod' x (2 ^ l); omega

/-- During the scan every key before idx has a high part at most h, and idx is inside or past the bucket. -/
def ScanInv (h idx : ℕ) : Prop := below (n := n) (l := l) (v := v) h ≤ idx ∧ ∀ j < idx, hi l v j ≤ h

/-- Inside the scan, bit h + idx is set exactly when key idx is in bucket h. -/
theorem H_bucket (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {h idx : ℕ} (hinv : ScanInv (n := n) (l := l) (v := v) h idx)
    (hidx : idx < n) : H n l v (h + idx) ↔ hi l v idx = h := by
  constructor
  · rintro ⟨j, hj, he⟩
    rcases lt_trichotomy j idx with hji | rfl | hji
    · have := hinv.2 j hji; omega
    · omega
    · have := (below_iff hs h hj).not.mpr (by have := hinv.1; omega); omega
  · intro he; exact ⟨idx, hidx, by omega⟩

/-- A scan only ever answers with a key equal to c. -/
theorem scan_sound (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) (h low : ℕ) :
    ∀ fuel idx i, ScanInv (n := n) (l := l) (v := v) h idx → scan (n := n) (l := l) (v := v) h low fuel idx = some i →
      i < n ∧ hi l v i = h ∧ lo l v i = low
  | 0, _, _, _, hsc => by simp [scan] at hsc
  | fuel + 1, idx, i, hinv, hsc => by
    unfold scan at hsc
    by_cases hc : idx < n ∧ H n l v (h + idx)
    · rw [if_pos hc] at hsc
      have hb := (H_bucket hs hinv hc.1).mp hc.2
      by_cases he : lo l v idx = low
      · rw [if_pos he] at hsc
        cases hsc
        exact ⟨hc.1, hb, he⟩
      · rw [if_neg he] at hsc
        by_cases hg : lo l v idx > low
        · rw [if_pos hg] at hsc; cases hsc
        · rw [if_neg hg] at hsc
          refine scan_sound hs h low fuel (idx + 1) i ⟨by have := hinv.1; omega, fun j hj => ?_⟩ hsc
          rcases Nat.lt_succ_iff_lt_or_eq.mp hj with hj | rfl
          · exact hinv.2 j hj
          · omega
    · rw [if_neg hc] at hsc; cases hsc

/-- A scan from the bucket start reaches key i when key i lies in the bucket and has the low part sought. -/
theorem scan_complete (hst : ∀ i j, i < j → j < n → v i < v j) (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j)
    {h i : ℕ} (hi' : i < n) (hh : hi l v i = h) :
    ∀ d fuel idx, idx + d = i → ScanInv (n := n) (l := l) (v := v) h idx → d < fuel →
      scan (n := n) (l := l) (v := v) h (lo l v i) fuel idx = some i
  | d, 0, _, _, _, hf => absurd hf (Nat.not_lt_zero _)
  | 0, fuel + 1, idx, he, hinv, _ => by
    simp only [Nat.add_zero] at he
    subst he
    have hH := (H_bucket hs hinv hi').mpr hh
    simp [scan, hi', hH]
  | d + 1, fuel + 1, idx, he, hinv, hf => by
    have hidx : idx < n := by omega
    have hmono := hi_mono (l := l) hs (show idx ≤ i by omega) hi'
    have hge : h ≤ hi l v idx := by
      by_contra hc
      have := (below_iff (l := l) hs h hidx).mp (by omega)
      have := hinv.1
      omega
    have hb : hi l v idx = h := by omega
    have hH := (H_bucket hs hinv hidx).mpr hb
    -- Same high part and a smaller key, so a smaller low part, and the scan moves on.
    have hlo : lo l v idx < lo l v i := by
      have hv := hst idx i (by omega) hi'
      have e1 := split_key (l := l) (v idx)
      have e2 := split_key (l := l) (v i)
      unfold hi at hb hh
      unfold lo
      rw [hb] at e1
      rw [hh] at e2
      omega
    unfold scan
    rw [if_pos ⟨hidx, hH⟩, if_neg (by omega), if_neg (by omega)]
    refine scan_complete hst hs hi' hh d fuel (idx + 1) (by omega) ⟨by have := hinv.1; omega, fun j hj => ?_⟩ (by omega)
    rcases Nat.lt_succ_iff_lt_or_eq.mp hj with hj | rfl
    · exact hinv.2 j hj
    · omega

/-- The scan starts at the bucket's first key. -/
theorem start_eq (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) {c p0 : ℕ}
    (hp0 : 1 ≤ c / 2 ^ l → Sel0 n l v (c / 2 ^ l - 1) p0) :
    (if c / 2 ^ l = 0 then 0 else p0 + 1) - c / 2 ^ l = below (n := n) (l := l) (v := v) (c / 2 ^ l) := by
  by_cases h0 : c / 2 ^ l = 0
  · rw [if_pos h0, h0]
    simp [below]
  · rw [if_neg h0]
    exact sel0_bucket (h := c / 2 ^ l) hs (Nat.one_le_iff_ne_zero.mpr h0) (hp0 (Nat.one_le_iff_ne_zero.mpr h0))

theorem start_inv (hs : ∀ i j, i ≤ j → j < n → v i ≤ v j) (h : ℕ) :
    ScanInv (n := n) (l := l) (v := v) h (below (n := n) (l := l) (v := v) h) :=
  ⟨le_rfl, fun j hj => le_of_lt ((below_iff hs h (by have := below_le (n := n) (l := l) (v := v) h; omega)).mpr hj)⟩

/-- `rank(c)` answers only with the index of c. -/
theorem rank_some (hst : ∀ i j, i < j → j < n → v i < v j) {c p0 i : ℕ}
    (hp0 : 1 ≤ c / 2 ^ l → Sel0 n l v (c / 2 ^ l - 1) p0) (h : rank (n := n) (l := l) (v := v) c p0 = some i) :
    i < n ∧ v i = c := by
  have hs : ∀ i j, i ≤ j → j < n → v i ≤ v j := fun i j hij hj => by
    rcases Nat.lt_or_eq_of_le hij with h | rfl
    · exact le_of_lt (hst i j h hj)
    · exact le_rfl
  unfold rank at h
  simp only at h
  rw [start_eq hs hp0] at h
  obtain ⟨hi', hh, hl⟩ := scan_sound hs _ _ _ _ i (start_inv hs _) h
  refine ⟨hi', ?_⟩
  rw [split_key (l := l) (v i), split_key (l := l) c]
  unfold hi at hh
  unfold lo at hl
  rw [hh, hl]

/-- `rank(c)` finds every stored key. -/
theorem rank_complete (hst : ∀ i j, i < j → j < n → v i < v j) {c p0 i : ℕ}
    (hp0 : 1 ≤ c / 2 ^ l → Sel0 n l v (c / 2 ^ l - 1) p0) (hi' : i < n) (hv : v i = c) :
    rank (n := n) (l := l) (v := v) c p0 = some i := by
  have hs : ∀ i j, i ≤ j → j < n → v i ≤ v j := fun i j hij hj => by
    rcases Nat.lt_or_eq_of_le hij with h | rfl
    · exact le_of_lt (hst i j h hj)
    · exact le_rfl
  unfold rank
  simp only
  rw [start_eq hs hp0]
  have hh : hi l v i = c / 2 ^ l := by unfold hi; rw [hv]
  have hl : lo l v i = c % 2 ^ l := by unfold lo; rw [hv]
  have hb : below (n := n) (l := l) (v := v) (c / 2 ^ l) ≤ i := by
    by_contra hc
    have := (below_iff (l := l) hs (c / 2 ^ l) hi').mpr (by omega)
    omega
  rw [← hl]
  exact scan_complete hst hs hi' hh (i - below (n := n) (l := l) (v := v) (c / 2 ^ l)) (n + 1) _ (by omega)
    (start_inv hs _) (by omega)

end OptikProofs.Ef
