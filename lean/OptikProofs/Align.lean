import Mathlib.Data.List.Basic
import Mathlib.Algebra.BigOperators.Group.List.Basic

/-!
# Left-aligning indels (Prop. 6 and 6′)

Models `left_align` in `tools/src/map/cigar.rs`.

An alignment is a list of columns.
An M column pairs a query base with a text base and costs `mm q t`, while an I or a D column costs 1 per base, which is the edit distance behind NM.
`mm` is left abstract, so the result holds for the code's rule that any non-ACGT base mismatches, and for any other per-column mismatch cost.

`left_align` moves a run of `len` deletions (or insertions) left past the j match columns before it.
`shift_deletion` and `shift_insertion` show that the shifted columns align the same two strings, that the gap costs the same, and that the total changes by exactly the code's `delta`.
That delta compares the j re-paired columns before and after the move.
The code accepts a shift only when delta ≤ 0, so a shift never worsens the alignment (Prop. 6).
The passed columns need not be exact matches, so the same statement covers a mismatch beside the gap that moves with it (Prop. 6′), as in the HG002 deletion-plus-SNP case.
-/

namespace OptikProofs.Align

variable {α : Type*} (mm : α → α → ℕ)

/-- CIGAR operations of an M/I/D alignment. -/
inductive Op
  | M
  | I
  | D
  deriving DecidableEq

/-- Total cost of the columns `ops` aligning `q` against `t`, reading both from the front. -/
def cost : List Op → List α → List α → ℕ
  | [], _, _ => 0
  | Op.M :: ops, a :: q, b :: t => mm a b + cost ops q t
  | Op.I :: ops, _ :: q, t => 1 + cost ops q t
  | Op.D :: ops, q, _ :: t => 1 + cost ops q t
  | _, _, _ => 0

/-- Total mismatch cost of pairing a query run with an equal-length text run, base by base. -/
def pairCost : List α → List α → ℕ
  | a :: q, b :: t => mm a b + pairCost q t
  | _, _ => 0

/-- A run of match columns pairs its bases one to one, then the remaining columns continue. -/
theorem cost_matches_then (rest : List Op) (q' t' : List α) :
    ∀ (a t : List α), a.length = t.length →
      cost mm (List.replicate a.length Op.M ++ rest) (a ++ q') (t ++ t') = pairCost mm a t + cost mm rest q' t'
  | [], [], _ => by simp [pairCost]
  | x :: a, y :: t, h => by
    simp only [List.length_cons, Nat.add_right_cancel_iff] at h
    simp only [List.length_cons, List.replicate_succ, List.cons_append, cost, pairCost,
      cost_matches_then rest q' t' a t h, Nat.add_assoc]

/-- A run of deletions consumes text only and costs one per base. -/
theorem cost_dels_then (rest : List Op) (q t' : List α) :
    ∀ (t : List α), cost mm (List.replicate t.length Op.D ++ rest) q (t ++ t') = t.length + cost mm rest q t'
  | [] => by simp
  | y :: t => by
    simp only [List.length_cons, List.replicate_succ, List.cons_append, cost, cost_dels_then rest q t' t]
    omega

/-- A run of insertions consumes query only and costs one per base. -/
theorem cost_ins_then (rest : List Op) (q' t : List α) :
    ∀ (q : List α), cost mm (List.replicate q.length Op.I ++ rest) (q ++ q') t = q.length + cost mm rest q' t
  | [] => by simp
  | x :: q => by
    simp only [List.length_cons, List.replicate_succ, List.cons_append, cost, cost_ins_then rest q' t q]
    omega

/-- Prop. 6 and 6′ for a deletion.
The j query bases `a` sit before a deletion of `len` text bases, against the text run `t₀` of j + len bases.
Before the shift `a` pairs with the first j bases of `t₀` and the last `len` are deleted.
After it, the first `len` are deleted and `a` pairs with the last j.
Both column lists consume exactly `a` and `t₀`, the deletion costs `len` either way, and the cost moves only by the change in `pairCost`, which is the code's delta. -/
theorem shift_deletion (a t₀ : List α) (len : ℕ) (hlen : a.length + len = t₀.length) :
    cost mm (List.replicate a.length Op.M ++ List.replicate len Op.D) a t₀ = pairCost mm a (t₀.take a.length) + len ∧
      cost mm (List.replicate len Op.D ++ List.replicate a.length Op.M) a t₀ = len + pairCost mm a (t₀.drop len) := by
  constructor
  · have hd : (t₀.drop a.length).length = len := by simp; omega
    have := cost_matches_then mm (List.replicate len Op.D) [] (t₀.drop a.length) a (t₀.take a.length)
      (by simp; omega)
    rw [List.append_nil, List.take_append_drop] at this
    rw [this]
    have h2 := cost_dels_then mm [] [] [] (t₀.drop a.length)
    simp only [List.append_nil, hd] at h2
    rw [h2]
    simp [cost]
  · have ht : (t₀.take len).length = len := by simp; omega
    have := cost_dels_then mm (List.replicate a.length Op.M) a (t₀.drop len) (t₀.take len)
    rw [ht, List.take_append_drop] at this
    rw [this]
    have h2 := cost_matches_then mm [] [] [] a (t₀.drop len) (by simp; omega)
    simp only [List.append_nil] at h2
    rw [h2]
    simp [cost]

/-- Prop. 6 and 6′ for an insertion, which mirrors the deletion with `len` query bases inserted next to j text bases `b`. -/
theorem shift_insertion (q₀ b : List α) (len : ℕ) (hlen : b.length + len = q₀.length) :
    cost mm (List.replicate b.length Op.M ++ List.replicate len Op.I) q₀ b = pairCost mm (q₀.take b.length) b + len ∧
      cost mm (List.replicate len Op.I ++ List.replicate b.length Op.M) q₀ b = len + pairCost mm (q₀.drop len) b := by
  constructor
  · have hd : (q₀.drop b.length).length = len := by simp; omega
    have hb : (q₀.take b.length).length = b.length := by simp; omega
    have := cost_matches_then mm (List.replicate len Op.I) (q₀.drop b.length) [] (q₀.take b.length) b hb
    rw [hb, List.take_append_drop, List.append_nil] at this
    rw [this]
    have h2 := cost_ins_then mm [] [] [] (q₀.drop b.length)
    simp only [List.append_nil, hd] at h2
    rw [h2]
    simp [cost]
  · have ht : (q₀.take len).length = len := by simp; omega
    have := cost_ins_then mm (List.replicate b.length Op.M) (q₀.drop len) b (q₀.take len)
    rw [ht, List.take_append_drop] at this
    rw [this]
    have hb : (q₀.drop len).length = b.length := by simp; omega
    have h2 := cost_matches_then mm [] [] [] (q₀.drop len) b hb
    rw [hb] at h2
    simp only [List.append_nil] at h2
    rw [h2]
    simp [cost]

/-- The code shifts only when the re-paired columns cost no more (delta ≤ 0), and then the alignment costs no more. -/
theorem shift_deletion_le (a t₀ : List α) (len : ℕ) (hlen : a.length + len = t₀.length)
    (hdelta : pairCost mm a (t₀.drop len) ≤ pairCost mm a (t₀.take a.length)) :
    cost mm (List.replicate len Op.D ++ List.replicate a.length Op.M) a t₀ ≤
      cost mm (List.replicate a.length Op.M ++ List.replicate len Op.D) a t₀ := by
  obtain ⟨h1, h2⟩ := shift_deletion mm a t₀ len hlen
  omega

theorem shift_insertion_le (q₀ b : List α) (len : ℕ) (hlen : b.length + len = q₀.length)
    (hdelta : pairCost mm (q₀.drop len) b ≤ pairCost mm (q₀.take b.length) b) :
    cost mm (List.replicate len Op.I ++ List.replicate b.length Op.M) q₀ b ≤
      cost mm (List.replicate b.length Op.M ++ List.replicate len Op.I) q₀ b := by
  obtain ⟨h1, h2⟩ := shift_insertion mm q₀ b len hlen
  omega

end OptikProofs.Align
