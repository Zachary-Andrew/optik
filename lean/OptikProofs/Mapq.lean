import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Algebra.Order.Round

/-!
# Short-read MAPQ as a posterior

The MAPQ code is `posterior_mapq_w` in `tools/src/map/mapq.rs`, with `edit_ratio` from `tools/src/map/paired.rs`.

Proposition 5 (prop:post) is the idealised identity.
Every candidate locus is equally likely beforehand, each edit multiplies a placement's likelihood by r, and the loci the seeds never reached add odds H relative to the best placement.
The best placement is then wrong with probability (S + H) / (1 + S + H), where S sums r^(dᵢ - d₀) over the competitors.

The code adds per-competitor weights, a floor on P(wrong), clamping, and an explicit MAPQ 0 rule.
Those are modelled separately as `implementedMapq`, so the paper can state both the identity and the formula the code runs.
The contrast theorem shows the MAPQ 0 for a tied competitor is that explicit rule, since the posterior alone gives MAPQ 3.
-/

namespace OptikProofs.Mapq

open Real

/-- P(wrong) given competitor odds `S` and unseen-locus odds `H`. -/
noncomputable def pWrong (S H : ℝ) : ℝ := (S + H) / (1 + S + H)

/-- Competitor odds S = Σ r^(dᵢ - d₀) for edit distances `ds` against the best placement's `d0`. -/
noncomputable def oddsSum (r : ℝ) (d0 : ℤ) (ds : List ℤ) : ℝ := (ds.map fun d => r ^ (d - d0)).sum

/-- prop:post.
With equal priors and likelihood c · r^d for a placement with d edits, the best placement's posterior of being wrong is `pWrong`.
The unseen loci contribute H times the best placement's likelihood. -/
theorem posterior_identity {r c H : ℝ} (hr : 0 < r) (hc : 0 < c) (hH : 0 ≤ H) (d0 : ℤ) (ds : List ℤ) :
    let L0 := c * r ^ d0
    let rest := (ds.map fun d => c * r ^ d).sum
    1 - L0 / (L0 + rest + H * L0) = pWrong (oddsSum r d0 ds) H := by
  intro L0 rest
  have hL0 : 0 < L0 := mul_pos hc (zpow_pos hr d0)
  have hrest : rest = L0 * oddsSum r d0 ds := by
    simp only [rest, L0, oddsSum, ← List.sum_map_mul_left]
    congr 1
    refine List.map_congr_left fun d _ => ?_
    rw [zpow_sub₀ hr.ne']
    field_simp
  have hS : 0 ≤ oddsSum r d0 ds := List.sum_nonneg fun x hx => by
    obtain ⟨d, -, rfl⟩ := List.mem_map.mp hx
    exact (zpow_pos hr _).le
  rw [hrest, pWrong]
  field_simp
  ring

theorem pWrong_nonneg {S H : ℝ} (hS : 0 ≤ S) (hH : 0 ≤ H) : 0 ≤ pWrong S H :=
  div_nonneg (add_nonneg hS hH) (by linarith)

theorem pWrong_lt_one {S H : ℝ} (hS : 0 ≤ S) (hH : 0 ≤ H) : pWrong S H < 1 := by
  rw [pWrong, div_lt_one (by linarith)]
  linarith

/-- More competitor odds can only make the best placement likelier to be wrong. -/
theorem pWrong_mono_S {S S' H : ℝ} (hS : 0 ≤ S) (hH : 0 ≤ H) (h : S ≤ S') : pWrong S H ≤ pWrong S' H := by
  unfold pWrong
  rw [div_le_div_iff₀ (by linarith) (by linarith)]
  nlinarith

/-- So can more unseen-locus odds. -/
theorem pWrong_mono_H {S H H' : ℝ} (hS : 0 ≤ S) (hH : 0 ≤ H) (h : H ≤ H') : pWrong S H ≤ pWrong S H' := by
  unfold pWrong
  rw [div_le_div_iff₀ (by linarith) (by linarith)]
  nlinarith

/-! ## From probability to MAPQ -/

/-- `SAM_MAX_MAPQ` in tools/src/map/sam.rs. -/
def samMaxMapq : ℤ := 60

/-- `P_WRONG_FLOOR` in tools/src/map/mapq.rs. -/
noncomputable def pWrongFloor : ℝ := 1e-9

/-- −10 log₁₀ p, rounded and clamped to [0, 60] as `posterior_mapq_w` does with `PHRED_SCALE = 10`. -/
noncomputable def phred (p : ℝ) : ℤ := max 0 (min samMaxMapq (round (-10 * logb 10 p)))

/-- A likelier error never earns a higher MAPQ. -/
theorem phred_antitone {p q : ℝ} (hp : 0 < p) (hpq : p ≤ q) : phred q ≤ phred p := by
  unfold phred
  have hlog : logb 10 p ≤ logb 10 q := logb_le_logb_of_le (by norm_num) hp hpq
  have hround : round (-10 * logb 10 q) ≤ round (-10 * logb 10 p) := by
    rw [round_eq, round_eq]
    exact Int.floor_le_floor (by linarith)
  omega

/-- MAPQ falls as competitor odds grow, for the idealised posterior. -/
theorem mapq_antitone_S {S S' H : ℝ} (hS : 0 ≤ S) (hH : 0 < H) (h : S ≤ S') :
    phred (pWrong S' H) ≤ phred (pWrong S H) :=
  phred_antitone (div_pos (by linarith) (by linarith)) (pWrong_mono_S hS hH.le h)

/-- MAPQ falls as unseen-locus odds grow. -/
theorem mapq_antitone_H {S H H' : ℝ} (hS : 0 ≤ S) (hH : 0 < H) (h : H ≤ H') :
    phred (pWrong S H') ≤ phred (pWrong S H) :=
  phred_antitone (div_pos (by linarith) (by linarith)) (pWrong_mono_H hS hH.le h)

/-! ## The formula as implemented -/

/-- `posterior_mapq_w` before its MAPQ 0 rule, flooring the weighted competitor odds plus the unseen-locus odds before the Phred conversion. -/
noncomputable def mapqNoOverride (best : ℤ) (comps : List (ℤ × ℝ)) (ratio hidden : ℝ) : ℤ :=
  let odds := (comps.map fun c => c.2 * ratio ^ (c.1 - best)).sum + hidden
  phred (max pWrongFloor (odds / (1 + odds)))

open Classical in
/-- `posterior_mapq_w`, `ratio` standing for `edit_ratio(best_ed, len)`.
A competitor at least as close as the best placement, with weight at least 1, sets MAPQ to 0 outright. -/
noncomputable def implementedMapq (best : ℤ) (comps : List (ℤ × ℝ)) (ratio hidden : ℝ) : ℤ :=
  if ∃ c ∈ comps, c.1 ≤ best ∧ 1 ≤ c.2 then 0 else mapqNoOverride best comps ratio hidden

theorem implementedMapq_override {best : ℤ} {comps : List (ℤ × ℝ)} {ratio hidden : ℝ}
    (h : ∃ c ∈ comps, c.1 ≤ best ∧ 1 ≤ c.2) : implementedMapq best comps ratio hidden = 0 := by
  rw [implementedMapq, if_pos h]

/-- −10 log₁₀ ½ rounds to 3, because 2¹⁰ ≥ 10³ and 2²⁰ < 10⁷ pin log₁₀ 2 inside [0.3, 0.35). -/
theorem phred_half : phred (1 / 2) = 3 := by
  have h2 : logb 10 (1 / 2) = -logb 10 2 := by
    rw [one_div, logb_inv]
  have hlo : (3 : ℝ) ≤ 10 * logb 10 2 := by
    have : logb 10 ((10 : ℝ) ^ 3) ≤ logb 10 ((2 : ℝ) ^ 10) :=
      logb_le_logb_of_le (by norm_num) (by norm_num) (by norm_num)
    rw [logb_pow, logb_pow, logb_self_eq_one (by norm_num)] at this
    push_cast at this
    linarith
  have hhi : 10 * logb 10 2 < (7 : ℝ) / 2 := by
    have : logb 10 ((2 : ℝ) ^ 20) < logb 10 ((10 : ℝ) ^ 7) :=
      logb_lt_logb (by norm_num) (by norm_num) (by norm_num)
    rw [logb_pow, logb_pow, logb_self_eq_one (by norm_num)] at this
    push_cast at this
    linarith
  have hround : round (-10 * logb 10 (1 / 2)) = 3 := by
    rw [h2, round_eq, Int.floor_eq_iff]
    constructor <;> push_cast <;> linarith
  rw [phred, hround]
  decide

/-- Without the explicit rule, one tied competitor of weight 1 and no unseen-locus odds give P(wrong) = ½ and MAPQ 3, not 0.
So the MAPQ 0 the paper reports for ties is the rule in the code, not a consequence of Proposition 5. -/
theorem tie_without_override (best : ℤ) (ratio : ℝ) :
    mapqNoOverride best [(best, 1)] ratio 0 = 3 := by
  have hmax : max pWrongFloor ((1 : ℝ) / (1 + 1)) = 1 / 2 := by
    rw [max_eq_right] <;> norm_num [pWrongFloor]
  simp only [mapqNoOverride, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, sub_self,
    zpow_zero, mul_one, add_zero]
  rw [hmax, phred_half]

/-- With the rule, the same tie gets MAPQ 0. -/
theorem tie_with_override (best : ℤ) (ratio hidden : ℝ) : implementedMapq best [(best, 1)] ratio hidden = 0 :=
  implementedMapq_override ⟨(best, 1), by simp, le_refl _, le_refl _⟩

/-! ## Calibration test (prop:cal, partial) -/

/-- The calibration test rests on Σ pᵢ² ≤ λ · max pᵢ, where λ = Σ pᵢ is the expected number of wrong placements in a MAPQ bin. -/
theorem sum_sq_le (ps : List ℝ) (m : ℝ) (hnn : ∀ p ∈ ps, 0 ≤ p) (hm : ∀ p ∈ ps, p ≤ m) :
    (ps.map (· ^ 2)).sum ≤ ps.sum * m := by
  induction ps with
  | nil => simp
  | cons p ps ih =>
    have hp : p ^ 2 ≤ p * m := by
      rw [sq]
      exact mul_le_mul_of_nonneg_left (hm p (by simp)) (hnn p (by simp))
    have hrest := ih (fun q hq => hnn q (by simp [hq])) (fun q hq => hm q (by simp [hq]))
    simp only [List.map_cons, List.sum_cons]
    linarith

/-- prop:cal with Le Cam's theorem as a named hypothesis.
If the distance between the bin's Poisson-binomial count and Poisson(λ) is at most Σ pᵢ², as Le Cam (1960) shows, it is also at most λ · max pᵢ.
Le Cam's theorem itself is cited rather than formalized. -/
theorem calibration_bound (ps : List ℝ) (m tv : ℝ) (hnn : ∀ p ∈ ps, 0 ≤ p) (hm : ∀ p ∈ ps, p ≤ m)
    (hLeCam : tv ≤ (ps.map (· ^ 2)).sum) : tv ≤ ps.sum * m :=
  hLeCam.trans (sum_sq_le ps m hnn hm)

end OptikProofs.Mapq
