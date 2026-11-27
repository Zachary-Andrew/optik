import Mathlib.Data.List.Sort
import Mathlib.Order.Interval.Finset.Nat
import OptikProofs.Kmer

/-!
# Window minimizers

Models `extract_minimizers` in `include/minimizer.hpp` and the m-mer minimizer that SpssCore picks inside a k-mer.

`h p` is the order value of the k-mer at position p, a hash in the code.
Window j covers positions j to j + w - 1, and its minimizer is the leftmost position holding the smallest value.
The streaming loop tracks that position with a strict `<` update and rescans only when the minimum slides out of the window, and `stream_isLeftMin` proves it always holds the leftmost minimum.
Leftmost ties are what make minimizer positions non-decreasing from one window to the next (`minPos_mono`), and the parallel seed build needs that order.

The model covers one run of valid bases.
The code restarts at every N, and since no window spans an N, each run matches this model separately.
-/

namespace OptikProofs.Minimizer

variable {α : Type*} [LinearOrder α] (h : ℕ → α)

/-- `i` is the leftmost position of [a, b] holding that range's smallest value. -/
def IsLeftMin (a b i : ℕ) : Prop :=
  a ≤ i ∧ i ≤ b ∧ (∀ j, a ≤ j → j ≤ b → h i ≤ h j) ∧ ∀ j, a ≤ j → j < i → h i < h j

theorem IsLeftMin.unique {a b i i' : ℕ} (hi : IsLeftMin h a b i) (hi' : IsLeftMin h a b i') : i = i' := by
  obtain ⟨ha, hb, hmin, hleft⟩ := hi
  obtain ⟨ha', hb', hmin', hleft'⟩ := hi'
  by_contra hne
  rcases Nat.lt_or_gt_of_ne hne with hlt | hlt
  · exact absurd (hmin i' ha' hb') (not_le.mpr (hleft' i ha hlt))
  · exact absurd (hmin' i ha hb) (not_le.mpr (hleft i' ha' hlt))

/-- The rescan in `extract_minimizers`, which walks [a, a + n] and moves to a position only when its value is strictly smaller. -/
def scanFrom (best : ℕ) : List ℕ → ℕ
  | [] => best
  | j :: js => scanFrom (if h j < h best then j else best) js

/-- The leftmost minimum of [a, b], computed the way the code rescans. -/
def scanMin (a b : ℕ) : ℕ := scanFrom h a (List.range' (a + 1) (b - a))

theorem scanFrom_spec (a : ℕ) :
    ∀ (n best : ℕ), IsLeftMin h a (a + n) best → ∀ m,
      IsLeftMin h a (a + n + m) (scanFrom h best (List.range' (a + n + 1) m))
  | n, best, hb, 0 => by simpa [scanFrom] using hb
  | n, best, hb, m + 1 => by
    rw [List.range'_succ, scanFrom]
    have key : IsLeftMin h a (a + (n + 1)) (if h (a + n + 1) < h best then a + n + 1 else best) := by
      obtain ⟨ha, hle, hmin, hleft⟩ := hb
      split_ifs with hlt
      · refine ⟨by omega, by omega, fun j hj hj' => ?_, fun j hj hj' => ?_⟩
        · rcases Nat.lt_or_ge j (a + n + 1) with hj'' | hj''
          · exact le_of_lt (lt_of_lt_of_le hlt (hmin j hj (by omega)))
          · rw [show j = a + n + 1 by omega]
        · exact lt_of_lt_of_le hlt (hmin j hj (by omega))
      · refine ⟨ha, by omega, fun j hj hj' => ?_, hleft⟩
        rcases Nat.lt_or_ge j (a + n + 1) with hj'' | hj''
        · exact hmin j hj (by omega)
        · rw [show j = a + n + 1 by omega]; exact not_lt.mp hlt
    have := scanFrom_spec a (n + 1) _ key m
    simpa [Nat.add_assoc, Nat.add_comm 1 m, Nat.add_left_comm] using this

theorem scanMin_spec {a b : ℕ} (hab : a ≤ b) : IsLeftMin h a b (scanMin h a b) := by
  have base : IsLeftMin h a (a + 0) a :=
    ⟨le_refl a, by omega, fun j hj hj' => by rw [show j = a by omega], fun j hj hj' => by omega⟩
  have := scanFrom_spec h a 0 a base (b - a)
  simpa [scanMin, Nat.add_sub_cancel' hab] using this

/-- The minimizer position of window j, which spans w ≥ 1 positions. -/
def minPos (w j : ℕ) : ℕ := scanMin h j (j + w - 1)

theorem minPos_spec {w : ℕ} (hw : 1 ≤ w) (j : ℕ) : IsLeftMin h j (j + w - 1) (minPos h w j) :=
  scanMin_spec h (by omega)

/-- Moving the window right never moves its leftmost minimum left.
A smaller position in the next window would also lie in this one with a value no larger, which leftmost ties rule out. -/
theorem minPos_mono {w : ℕ} (hw : 1 ≤ w) {j : ℕ} : minPos h w j ≤ minPos h w (j + 1) := by
  obtain ⟨ha, hb, hmin, hleft⟩ := minPos_spec h hw j
  obtain ⟨ha', hb', hmin', -⟩ := minPos_spec h hw (j + 1)
  by_contra hlt
  replace hlt := not_le.mp hlt
  have h1 : h (minPos h w (j + 1)) ≤ h (minPos h w j) := hmin' _ (by omega) (by omega)
  have h2 : h (minPos h w j) < h (minPos h w (j + 1)) := hleft _ (by omega) hlt
  exact absurd h1 (not_le.mpr h2)

theorem minPos_monotone {w : ℕ} (hw : 1 ≤ w) : Monotone (minPos h w) :=
  monotone_nat_of_le_succ fun _ => minPos_mono h hw

/-! ## The streaming loop -/

/-- In `extract_minimizers` a new position replaces the tracked one only with a strictly smaller value, so a tie leaves the earlier position tracked. -/
def update (idx cur : ℕ) : ℕ := if cur = 0 then 0 else if h cur < h idx then cur else idx

/-- One step at position `cur`, given the tracked position `idx`.
A minimum that slid out of the window triggers a rescan. -/
def step (w idx cur : ℕ) : ℕ :=
  if w ≤ cur + 1 ∧ update h idx cur < cur + 1 - w then scanMin h (cur + 1 - w) cur else update h idx cur

/-- The tracked position after positions 0 to `cur` have been read. -/
def stream (w : ℕ) : ℕ → ℕ
  | 0 => step h w 0 0
  | cur + 1 => step h w (stream w cur) (cur + 1)

/-- After position `cur` the loop holds the leftmost minimum of the positions read so far, limited to the current window once it is full. -/
theorem stream_isLeftMin {w : ℕ} (hw : 1 ≤ w) :
    ∀ cur, IsLeftMin h (cur + 1 - w) cur (stream h w cur) := by
  intro cur
  induction cur with
  | zero =>
    have hu : update h 0 0 = 0 := by simp [update]
    have hs : stream h w 0 = 0 := by
      simp only [stream, step, hu]
      rw [if_neg (by omega)]
    rw [hs, show (0 : ℕ) + 1 - w = 0 by omega]
    exact ⟨le_refl 0, le_refl 0, fun j hj hj' => by rw [show j = 0 by omega], fun j hj hj' => by omega⟩
  | succ cur ih =>
    obtain ⟨ha, hb, hmin, hleft⟩ := ih
    set idx := stream h w cur
    -- The candidate after the strict update is the leftmost minimum of [cur + 1 - w, cur + 1].
    set c := if h (cur + 1) < h idx then cur + 1 else idx with hc_def
    have upd : IsLeftMin h (cur + 1 - w) (cur + 1) c := by
      rw [hc_def]
      split_ifs with hlt
      · refine ⟨by omega, le_refl _, fun j hj hj' => ?_, fun j hj hj' => ?_⟩
        · rcases Nat.lt_or_ge j (cur + 1) with hj'' | hj''
          · exact le_of_lt (lt_of_lt_of_le hlt (hmin j hj (by omega)))
          · rw [show j = cur + 1 by omega]
        · exact lt_of_lt_of_le hlt (hmin j hj (by omega))
      · refine ⟨ha, by omega, fun j hj hj' => ?_, hleft⟩
        rcases Nat.lt_or_ge j (cur + 1) with hj'' | hj''
        · exact hmin j hj (by omega)
        · rw [show j = cur + 1 by omega]; exact not_lt.mp hlt
    have hu : update h idx (cur + 1) = c := by simp [update, hc_def]
    have hstep : step h w idx (cur + 1) =
        if w ≤ cur + 1 + 1 ∧ c < cur + 1 + 1 - w then scanMin h (cur + 1 + 1 - w) (cur + 1) else c := by
      simp only [step, hu]
    show IsLeftMin h (cur + 1 + 1 - w) (cur + 1) (step h w idx (cur + 1))
    rw [hstep]
    split_ifs with hslide
    · exact scanMin_spec h (by omega)
    · -- Without a rescan the candidate still lies in the window, so it is still leftmost in the narrower range.
      obtain ⟨ua, ub, umin, uleft⟩ := upd
      have hc : cur + 1 + 1 - w ≤ c := by
        by_contra hc
        exact hslide ⟨by omega, by omega⟩
      exact ⟨hc, ub, fun j hj hj' => umin j (by omega) hj', fun j hj hj' => uleft j (by omega) hj'⟩

/-- Once a window is full, the streaming loop reports exactly that window's minimizer. -/
theorem stream_eq_minPos {w : ℕ} (hw : 1 ≤ w) (j : ℕ) : stream h w (j + w - 1) = minPos h w j := by
  have hs := stream_isLeftMin h hw (j + w - 1)
  rw [show j + w - 1 + 1 - w = j by omega] at hs
  exact hs.unique h (minPos_spec h hw j)

/-! ## Mirrored ties under reverse complement (Lemma 4.2) -/

/-- Offsets in [0, L] holding the smallest value, all of which a lookup tries. -/
def Tied (f : ℕ → α) (L o : ℕ) : Prop := o ≤ L ∧ ∀ o', o' ≤ L → f o ≤ f o'

/-- Reading offsets back to front mirrors the tied set, so o is tied for `f (L - ·)` exactly when L - o is tied for `f`. -/
theorem tied_mirror (f : ℕ → α) {L o : ℕ} (ho : o ≤ L) : Tied (fun o' => f (L - o')) L o ↔ Tied f L (L - o) := by
  constructor
  · rintro ⟨-, hmin⟩
    refine ⟨by omega, fun o' ho' => ?_⟩
    have := hmin (L - o') (by omega)
    simp only at this
    rwa [show L - (L - o') = o' by omega] at this
  · rintro ⟨-, hmin⟩
    exact ⟨ho, fun o' ho' => hmin (L - o') (by omega)⟩

open OptikProofs.Kmer in
/-- The m-mer at offset o of rc x is the reverse complement of x's m-mer at offset L - o, where L = k - m. -/
theorem mmer_rc {x : List Base} {m o : ℕ} (ho : o + m ≤ x.length) :
    ((rc x).drop o).take m = rc ((x.drop (x.length - m - o)).take m) := by
  apply List.ext_getElem
  · simp [rc]; omega
  · intro i h1 h2
    simp only [List.getElem_take, List.getElem_drop, rc, List.getElem_reverse, List.getElem_map]
    congr 2
    simp at h1 h2 ⊢
    omega

omit [LinearOrder α] in
open OptikProofs.Kmer in
/-- Lemma 4.2.
When the order depends only on the canonical m-mer, the order values of rc x read back to front are those of x.
So the tied minimizer offsets of rc x are the mirror image of x's, and a lookup trying every tie finds the same records on either strand. -/
theorem order_mirror (ord : ℕ → α) (x : List Base) {m : ℕ} (hm : m ≤ x.length) {o : ℕ}
    (ho : o ≤ x.length - m) :
    ord (canon (((rc x).drop o).take m)) = ord (canon ((x.drop (x.length - m - o)).take m)) := by
  rw [mmer_rc (by omega), canon_rc]

end OptikProofs.Minimizer
