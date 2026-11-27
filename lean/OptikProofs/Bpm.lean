import Std.Tactic.BVDecide

/-!
# Bit-parallel semi-global edit distance equals the DP

`bpm_semiglobal_ed` in `include/bpm_edit.hpp` is the verifier behind every short-read placement.

The spec is the unit-cost DP with the query end to end and the text ends free.
`D j r` is the cost of the first r query bases against the best text substring ending after text base j, so `D j 0 = 0` (free start) and `D 0 r = r`.
The answer is the least `D j m` over the columns j.

The code stores each column as vertical deltas, one bit per query row in the words `Pv` and `Mv`, and advances 64 rows per word in a few operations (Myers 1999, Hyyrö 2003 for several words).
Rows past the query's end in the last word compare as mismatches, which is what their all-zero match bits give, and no row below them depends on them.

* `cell` checks the per-row delta rule against the DP's minimum in all 32 cases.
* `carry` (by `bv_decide`) shows that the addition in the word step is the chain the rule needs, each row's negative horizontal delta feeding the next row.
* `column` puts them together over all rows and words, by induction on the row.
* `bpm_eq_dp` concludes that the code returns the least `D j m`, capped at `max_ed + 1`.
-/

namespace OptikProofs.Bpm

/-! ## The DP -/

variable (cost : Nat → Nat → Int)

/-- The DP over text columns j and query rows r, where `cost j r` is 0 when text base j matches query row r and 1 otherwise. -/
def D : Nat → Nat → Int
  | 0, r => r
  | _ + 1, 0 => 0
  | j + 1, r + 1 => min (min (D j r + cost j r) (D j (r + 1) + 1)) (D (j + 1) r + 1)
termination_by j r => (j, r)

/-- The vertical delta into row r + 1 of column j. -/
def vd (j r : Nat) : Int := D cost j (r + 1) - D cost j r

/-- The horizontal delta at row r from column j to column j + 1. -/
def hd (j r : Nat) : Int := D cost (j + 1) r - D cost j r

variable {cost}

/-- One DP cell in delta form, where the new cell is the old corner plus the least of the three moves. -/
theorem row (j r : Nat) :
    hd cost j (r + 1) = min (cost j r) (min (vd cost j r + 1) (hd cost j r + 1)) - vd cost j r ∧
      vd cost (j + 1) r = min (cost j r) (min (vd cost j r + 1) (hd cost j r + 1)) - hd cost j r := by
  simp only [hd, vd]
  rw [D]
  constructor <;> omega

theorem hd_zero (j : Nat) : hd cost j 0 = 0 := by cases j <;> simp [hd, D]

/-- Every delta is -1, 0 or 1. -/
theorem bounds (hc : ∀ j r, cost j r = 0 ∨ cost j r = 1) :
    ∀ j r, (-1 ≤ vd cost j r ∧ vd cost j r ≤ 1) ∧ (-1 ≤ hd cost j r ∧ hd cost j r ≤ 1) := by
  intro j
  induction j with
  | zero =>
    intro r
    induction r with
    | zero => simp [vd, hd, D]
    | succ r ih =>
      refine ⟨by simp [vd, D]; omega, ?_⟩
      have := (row (cost := cost) 0 r).1
      have := hc 0 r
      have hv : vd cost 0 r = 1 := by simp [vd, D]; omega
      omega
  | succ j ihj =>
    intro r
    induction r with
    | zero =>
      have h := (row (cost := cost) j 0).2
      have := hc j 0
      have := (ihj 0).1
      rw [hd_zero] at h
      refine ⟨by omega, by rw [hd_zero]; omega⟩
    | succ r ih =>
      have h1 := (row (cost := cost) j (r + 1)).2
      have h2 := (row (cost := cost) (j + 1) r).1
      have := hc j (r + 1)
      have := hc (j + 1) r
      have := (ihj (r + 1)).1
      have := (ihj (r + 1)).2
      have := ih.1
      have := ih.2
      have hrow := (row (cost := cost) j r).1
      have := (ihj r).1
      have := (ihj r).2
      refine ⟨by omega, by omega⟩

/-! ## One row of bits -/

def b2i (b : Bool) : Int := if b then 1 else 0

/-- The word step's formulas, one row at a time, give the DP's new deltas.
`eq` is the match bit, (pv, mv) the old vertical delta and (phin, mhin) the horizontal delta of the row below. -/
theorem cell : ∀ eq pv mv phin mhin : Bool, (pv && mv) = false → (phin && mhin) = false →
    let xh := eq || mhin
    let ph := mv || !(xh || pv)
    let mh := pv && xh
    let xv := eq || mv
    let pv' := mhin || !(xv || phin)
    let mv' := phin && xv
    let mn := min (if eq then 0 else 1) (min (b2i pv - b2i mv + 1) (b2i phin - b2i mhin + 1))
    b2i ph - b2i mh = mn - (b2i pv - b2i mv) ∧ b2i pv' - b2i mv' = mn - (b2i phin - b2i mhin) ∧
      (ph && mh) = false ∧ (pv' && mv') = false := by
  decide

/-! ## The word step -/

/-- The words one text base produces from one query word, transcribing the loop body of `bpm_semiglobal_ed`. -/
structure Out where
  pv : BitVec 64
  mv : BitVec 64
  ph : BitVec 64
  mh : BitVec 64

def wordStep (E P M : BitVec 64) (hin : Int) : Out :=
  let hneg : BitVec 64 := if hin < 0 then 1 else 0
  let hpos : BitVec 64 := if hin > 0 then 1 else 0
  let xv := E ||| M
  let E' := E ||| hneg
  let xh := (((E' &&& P) + P) ^^^ P) ||| E'
  let ph := M ||| ~~~(xh ||| P)
  let mh := P &&& xh
  let phs := (ph <<< 1) ||| hpos
  let mhs := (mh <<< 1) ||| hneg
  ⟨mhs ||| ~~~(xv ||| phs), phs &&& xv, ph, mh⟩

/-- The addition solves the chain in which each row's match run carries into the next row. -/
theorem carry (E P : BitVec 64) :
    (((E &&& P) + P) ^^^ P) ||| E = E ||| ((P &&& ((((E &&& P) + P) ^^^ P) ||| E)) <<< 1) := by
  bv_decide

theorem getLsbD_ite (c : Prop) [Decidable c] (b : Nat) :
    (if c then (1 : BitVec 64) else 0).getLsbD b = (decide c && decide (b = 0)) := by
  split <;> simp_all [BitVec.getLsbD_one]

/-- Row b of the word step reads the bit under it, or at bit 0 the delta coming in from the word before. -/
def below (x : BitVec 64) (top : Bool) (b : Nat) : Bool := if b = 0 then top else x.getLsbD (b - 1)

/-- Each output bit of the word step is the `cell` formula applied to its row. -/
theorem wordStep_bits (E P M : BitVec 64) (hin : Int) {b : Nat} (hb : b < 64) :
    let o := wordStep E P M hin
    let phin := below o.ph (decide (hin > 0)) b
    let mhin := below o.mh (decide (hin < 0)) b
    o.ph.getLsbD b = (M.getLsbD b || !((E.getLsbD b || mhin) || P.getLsbD b)) ∧
      o.mh.getLsbD b = (P.getLsbD b && (E.getLsbD b || mhin)) ∧
      o.pv.getLsbD b = (mhin || !((E.getLsbD b || M.getLsbD b) || phin)) ∧
      o.mv.getLsbD b = (phin && (E.getLsbD b || M.getLsbD b)) := by
  intro o phin mhin
  have hx := carry (E ||| if hin < 0 then 1 else 0) P
  generalize hX : (((((E ||| if hin < 0 then 1 else 0) &&& P) + P) ^^^ P) ||| (E ||| if hin < 0 then 1 else 0)) = xh at hx
  have hmh : o.mh = P &&& xh := by rw [← hX]; rfl
  have hph : o.ph = M ||| ~~~(xh ||| P) := by rw [← hX]; rfl
  -- Row b's xh is its match bit or the negative delta of the row below.
  have hxb : xh.getLsbD b = (E.getLsbD b || mhin) := by
    rw [hx]
    simp only [mhin, below, hmh, BitVec.getLsbD_or, BitVec.getLsbD_shiftLeft, BitVec.getLsbD_and, getLsbD_ite]
    by_cases h0 : b = 0 <;> simp [h0, hb] <;> omega
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [hph, BitVec.getLsbD_or, BitVec.getLsbD_not, BitVec.getLsbD_or, hxb]; simp [hb]
  · rw [hmh, BitVec.getLsbD_and, hxb]
  · have hpv : o.pv = ((o.mh <<< 1) ||| (if hin < 0 then (1 : BitVec 64) else 0)) |||
        ~~~((E ||| M) ||| ((o.ph <<< 1) ||| (if hin > 0 then (1 : BitVec 64) else 0))) := rfl
    rw [hpv]
    simp only [phin, mhin, below, BitVec.getLsbD_or, BitVec.getLsbD_not, BitVec.getLsbD_shiftLeft, getLsbD_ite]
    by_cases h0 : b = 0 <;> simp [h0, hb] <;> omega
  · have hmv : o.mv = ((o.ph <<< 1) ||| (if hin > 0 then (1 : BitVec 64) else 0)) &&& (E ||| M) := rfl
    rw [hmv]
    simp only [phin, below, BitVec.getLsbD_or, BitVec.getLsbD_and, BitVec.getLsbD_shiftLeft, getLsbD_ite]
    by_cases h0 : b = 0 <;> simp [h0, hb] <;> omega

/-! ## A column over several words -/

/-- The delta leaving a word at bit `bit`, which the code calls `hout`. -/
def hout (o : Out) (bit : Nat) : Int := if o.ph.getLsbD bit then 1 else if o.mh.getLsbD bit then -1 else 0

section Column

variable (E P M : Nat → BitVec 64)

/-- The horizontal delta entering word w is 0 for the first word (free text start) and otherwise what left bit 63 of the word before. -/
def hinOf : Nat → Int
  | 0 => 0
  | w + 1 => hout (wordStep (E w) (P w) (M w) (hinOf w)) 63

/-- Word w of the new column. -/
def outOf (w : Nat) : Out := wordStep (E w) (P w) (M w) (hinOf E P M w)

end Column

/-- Words `A` and `B` hold the deltas `d` of rows 0 .. R - 1, row r at bit r % 64 of word r / 64. -/
def Enc (A B : Nat → BitVec 64) (d : Nat → Int) (R : Nat) : Prop :=
  ∀ r, r < R → b2i ((A (r / 64)).getLsbD (r % 64)) - b2i ((B (r / 64)).getLsbD (r % 64)) = d r ∧
    ((A (r / 64)).getLsbD (r % 64) && (B (r / 64)).getLsbD (r % 64)) = false

theorem hout_spec {o : Out} {bit : Nat} (h : (o.ph.getLsbD bit && o.mh.getLsbD bit) = false) :
    decide (hout o bit > 0) = o.ph.getLsbD bit ∧ decide (hout o bit < 0) = o.mh.getLsbD bit ∧
      hout o bit = b2i (o.ph.getLsbD bit) - b2i (o.mh.getLsbD bit) := by
  unfold hout b2i
  cases h1 : o.ph.getLsbD bit <;> cases h2 : o.mh.getLsbD bit <;> simp_all

theorem b2i_decide (x : Int) (h : -1 ≤ x ∧ x ≤ 1) : b2i (decide (x > 0)) - b2i (decide (x < 0)) = x := by
  unfold b2i; by_cases h1 : x > 0 <;> by_cases h2 : x < 0 <;> simp [h1, h2] <;> omega

/-- Words holding column j's vertical deltas step to words holding column j + 1's, and each word's horizontal outputs are the DP's horizontal deltas. -/
theorem column (hc : ∀ j r, cost j r = 0 ∨ cost j r = 1) {E P M : Nat → BitVec 64} {j W : Nat}
    (hE : ∀ r, r < 64 * W → (E (r / 64)).getLsbD (r % 64) = decide (cost j r = 0))
    (henc : Enc P M (vd cost j) (64 * W)) :
    Enc (fun w => (outOf E P M w).ph) (fun w => (outOf E P M w).mh) (fun r => hd cost j (r + 1)) (64 * W) ∧
      Enc (fun w => (outOf E P M w).pv) (fun w => (outOf E P M w).mv) (vd cost (j + 1)) (64 * W) := by
  have key : ∀ r, r < 64 * W →
      (b2i ((outOf E P M (r / 64)).ph.getLsbD (r % 64)) - b2i ((outOf E P M (r / 64)).mh.getLsbD (r % 64)) = hd cost j (r + 1) ∧
        ((outOf E P M (r / 64)).ph.getLsbD (r % 64) && (outOf E P M (r / 64)).mh.getLsbD (r % 64)) = false) ∧
      (b2i ((outOf E P M (r / 64)).pv.getLsbD (r % 64)) - b2i ((outOf E P M (r / 64)).mv.getLsbD (r % 64)) = vd cost (j + 1) r ∧
        ((outOf E P M (r / 64)).pv.getLsbD (r % 64) && (outOf E P M (r / 64)).mv.getLsbD (r % 64)) = false) := by
    intro r
    induction r using Nat.strongRecOn with
    | ind r ih =>
      intro hr
      have hb : r % 64 < 64 := Nat.mod_lt _ (by decide)
      obtain ⟨hph, hmh, hpv, hmv⟩ := wordStep_bits (E (r / 64)) (P (r / 64)) (M (r / 64)) (hinOf E P M (r / 64)) hb
      -- Below the first row of a word is the incoming delta, and below any other row is the bit under it.
      have hin : b2i (below (outOf E P M (r / 64)).ph (decide (hinOf E P M (r / 64) > 0)) (r % 64)) -
            b2i (below (outOf E P M (r / 64)).mh (decide (hinOf E P M (r / 64) < 0)) (r % 64)) = hd cost j r ∧
          (below (outOf E P M (r / 64)).ph (decide (hinOf E P M (r / 64) > 0)) (r % 64) &&
            below (outOf E P M (r / 64)).mh (decide (hinOf E P M (r / 64) < 0)) (r % 64)) = false := by
        unfold below
        by_cases h0 : r % 64 = 0
        · simp only [h0, if_true]
          cases hw : r / 64 with
          | zero =>
            have : r = 0 := by omega
            subst this
            simp [hinOf, b2i, hd_zero]
          | succ w =>
            obtain ⟨⟨e1, e2⟩, -⟩ := ih (r - 1) (by omega) (by omega)
            rw [show (r - 1) / 64 = w by omega, show (r - 1) % 64 = 63 by omega, show r - 1 + 1 = r by omega] at e1
            rw [show (r - 1) / 64 = w by omega, show (r - 1) % 64 = 63 by omega] at e2
            obtain ⟨d1, d2, d3⟩ := hout_spec (o := outOf E P M w) (bit := 63) e2
            have hh : hinOf E P M (w + 1) = hout (outOf E P M w) 63 := rfl
            rw [hh, d1, d2]
            exact ⟨e1, e2⟩
        · simp only [h0, if_false]
          obtain ⟨⟨e1, e2⟩, -⟩ := ih (r - 1) (by omega) (by omega)
          rw [show (r - 1) / 64 = r / 64 by omega, show (r - 1) % 64 = r % 64 - 1 by omega,
            show r - 1 + 1 = r by omega] at e1
          rw [show (r - 1) / 64 = r / 64 by omega, show (r - 1) % 64 = r % 64 - 1 by omega] at e2
          exact ⟨e1, e2⟩
      obtain ⟨hv, hvv⟩ := henc r hr
      have hrow := row (cost := cost) j r
      have hcr := hc j r
      have hcost : (if decide (cost j r = 0) = true then (0 : Int) else 1) = cost j r := by
        rcases hcr with h | h <;> simp [h]
      have hcell := cell (E (r / 64) |>.getLsbD (r % 64)) ((P (r / 64)).getLsbD (r % 64)) ((M (r / 64)).getLsbD (r % 64))
        _ _ hvv hin.2
      rw [hE r hr, hcost, hv, hin.1] at hcell
      obtain ⟨c1, c2, c3, c4⟩ := hcell
      have eo : wordStep (E (r / 64)) (P (r / 64)) (M (r / 64)) (hinOf E P M (r / 64)) = outOf E P M (r / 64) := rfl
      rw [eo] at hph hmh hpv hmv
      rw [hE r hr] at hph hmh hpv hmv
      refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩⟩
      · rw [hph, hmh, c1, ← hrow.1]
      · rw [hph, hmh]; exact c3
      · rw [hpv, hmv, c2, ← hrow.2]
      · rw [hpv, hmv]; exact c4
  exact ⟨fun r hr => (key r hr).1, fun r hr => (key r hr).2⟩

/-! ## The text loop -/

/-- After some text bases the loop holds the delta words, the bottom-row score and the best score so far. -/
structure St where
  P : Nat → BitVec 64
  M : Nat → BitVec 64
  score : Int
  best : Int

/-- A text base with W words for an m-row query moves the bottom-row score by the delta leaving row m - 1. -/
def step (E : Nat → BitVec 64) (W m : Nat) (st : St) : St :=
  let score := st.score + hout (outOf E st.P st.M (W - 1)) ((m - 1) % 64)
  ⟨fun w => (outOf E st.P st.M w).pv, fun w => (outOf E st.P st.M w).mv, score, if score < st.best then score else st.best⟩

/-- `bpm_semiglobal_ed` after n text bases, with `E j` the match words of text base j. -/
def run (E : Nat → Nat → BitVec 64) (W m : Nat) : Nat → St
  | 0 => ⟨fun _ => BitVec.allOnes 64, fun _ => 0, m, m⟩
  | j + 1 => step (E j) W m (run E W m j)

variable (cost) in
/-- The spec's answer after n text bases. -/
def bestD (m : Nat) : Nat → Int
  | 0 => D cost 0 m
  | j + 1 => min (bestD m j) (D cost (j + 1) m)

/-- Myers' algorithm computes the DP.
After n text bases the score is `D n m`, the best is the least `D j m`, and the words hold column n. -/
theorem run_eq_dp (hc : ∀ j r, cost j r = 0 ∨ cost j r = 1) {E : Nat → Nat → BitVec 64} {W m : Nat}
    (hm : 0 < m) (hW : m ≤ 64 * W) (hW' : 64 * (W - 1) < m)
    (hE : ∀ j r, r < 64 * W → (E j (r / 64)).getLsbD (r % 64) = decide (cost j r = 0)) :
    ∀ n, (run E W m n).score = D cost n m ∧ (run E W m n).best = bestD cost m n ∧
      Enc (run E W m n).P (run E W m n).M (vd cost n) (64 * W)
  | 0 => by
    refine ⟨by simp [run, D], by simp [run, bestD, D], fun r _ => ?_⟩
    have : r % 64 < 64 := Nat.mod_lt _ (by decide)
    show b2i ((BitVec.allOnes 64).getLsbD (r % 64)) - b2i ((0 : BitVec 64).getLsbD (r % 64)) = vd cost 0 r ∧
      ((BitVec.allOnes 64).getLsbD (r % 64) && (0 : BitVec 64).getLsbD (r % 64)) = false
    rw [BitVec.getLsbD_allOnes]
    simp [b2i, vd, D, this]
    omega
  | n + 1 => by
    obtain ⟨hs, hbst, henc⟩ := run_eq_dp hc hm hW hW' hE n
    obtain ⟨hh, hv⟩ := column hc (hE n) henc
    obtain ⟨e1, e2⟩ := hh (m - 1) (by omega)
    simp only at e1 e2
    rw [show (m - 1) / 64 = W - 1 by omega, show m - 1 + 1 = m by omega] at e1
    rw [show (m - 1) / 64 = W - 1 by omega] at e2
    have hsc : (run E W m (n + 1)).score = D cost (n + 1) m := by
      simp only [run, step]
      rw [hs, (hout_spec e2).2.2, e1, hd]
      omega
    refine ⟨hsc, ?_, hv⟩
    have hb : (run E W m (n + 1)).best =
        if (run E W m (n + 1)).score < (run E W m n).best then (run E W m (n + 1)).score else (run E W m n).best := rfl
    rw [hb, hsc, hbst, bestD]
    split <;> omega

/-- The code's answer is the least `D j m` over the text, capped at `max_ed + 1`. -/
theorem bpm_eq_dp (hc : ∀ j r, cost j r = 0 ∨ cost j r = 1) {E : Nat → Nat → BitVec 64} {W m n : Nat} (max_ed : Int)
    (hm : 0 < m) (hW : m ≤ 64 * W) (hW' : 64 * (W - 1) < m)
    (hE : ∀ j r, r < 64 * W → (E j (r / 64)).getLsbD (r % 64) = decide (cost j r = 0)) :
    (if (run E W m n).best > max_ed then max_ed + 1 else (run E W m n).best) = min (bestD cost m n) (max_ed + 1) := by
  rw [(run_eq_dp hc hm hW hW' hE n).2.1]
  split <;> omega

end OptikProofs.Bpm
