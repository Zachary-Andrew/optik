import Mathlib.Order.MinMax
import Mathlib.Tactic.IntervalCases
import Mathlib.Data.List.Induction
import Std.Tactic.BVDecide

/-!
# K-mers, reverse complement and canonical form

The definitions follow `include/kmer.hpp`.
Bases are 2-bit codes with A = 0, C = 1, G = 2 and T = 3, so complementary bases are bitwise opposites.
A k-mer packs its first base into the highest two bits, which makes numeric order equal lexicographic order.

The list model states what the code computes.
The `BitVec` section then proves, for every k up to 32, that the word-level trick in `reverse_complement` computes that same reverse complement base by base.
-/

namespace OptikProofs.Kmer

/-- A 2-bit base code. -/
abbrev Base := Fin 4

/-- `base_complement`, pairing A with T and C with G. -/
def comp (b : Base) : Base := ⟨3 - b.val, by omega⟩

/-- The code complements with `base ^ kComplementMask`, and XOR with 3 is the same map on 2-bit codes. -/
theorem comp_val_eq_xor : ∀ b : Base, (comp b).val = b.val ^^^ 3 := by decide

@[simp] theorem comp_comp (b : Base) : comp (comp b) = b := by
  revert b; decide

/-- Reverse complement of a k-mer read 5' to 3'. -/
def rc (x : List Base) : List Base := (x.map comp).reverse

@[simp] theorem length_rc (x : List Base) : (rc x).length = x.length := by
  simp [rc]

/-- Reading the other strand twice gives the original k-mer back. -/
@[simp] theorem rc_rc (x : List Base) : rc (rc x) = x := by
  simp [rc, List.map_reverse, List.map_map, Function.comp_def]

/-- `encode_kmer`, which puts the first base in the highest two bits. -/
def pack (x : List Base) : ℕ := x.foldl (fun acc b => acc * 4 + b.val) 0

theorem pack_append (x : List Base) (b : Base) : pack (x ++ [b]) = pack x * 4 + b.val := by
  simp [pack, List.foldl_append]

theorem pack_lt (x : List Base) : pack x < 4 ^ x.length := by
  induction x using List.reverseRecOn with
  | nil => simp [pack]
  | append_singleton x b ih =>
    rw [pack_append, List.length_append, List.length_singleton, pow_succ]
    have := b.isLt
    omega

/-- Two k-mers of one length with the same packed word are the same k-mer, so 2k bits lose nothing. -/
theorem pack_injective {x y : List Base} (hlen : x.length = y.length) (h : pack x = pack y) : x = y := by
  induction x using List.reverseRecOn generalizing y with
  | nil =>
    simp at hlen
    exact (List.length_eq_zero_iff.mp hlen.symm).symm
  | append_singleton x b ih =>
    obtain ⟨y, c, rfl⟩ : ∃ y' c, y = y' ++ [c] := by
      rcases List.eq_nil_or_concat y with hy | ⟨y', c, hy⟩
      · subst hy; simp at hlen
      · exact ⟨y', c, by simpa using hy⟩
    simp only [List.length_append, List.length_singleton, Nat.add_right_cancel_iff] at hlen
    rw [pack_append, pack_append] at h
    have hb := b.isLt
    have hc := c.isLt
    have hbc : b.val = c.val := by omega
    have hxy : pack x = pack y := by omega
    rw [ih hlen hxy, Fin.ext hbc]

/-- The canonical form is the strand with the smaller packed word, as `canonical` does. -/
def canon (x : List Base) : ℕ := min (pack x) (pack (rc x))

/-- A k-mer and its reverse complement share one canonical form, so lookups are strand-blind. -/
@[simp] theorem canon_rc (x : List Base) : canon (rc x) = canon x := by
  simp [canon, min_comm]

/-- Equal canonical forms of same-length k-mers mean the same k-mer up to strand, and nothing else collides. -/
theorem canon_eq_iff {x y : List Base} (hlen : x.length = y.length) :
    canon x = canon y ↔ y = x ∨ y = rc x := by
  constructor
  · intro h
    have hrc : (rc x).length = y.length := by simp [hlen]
    simp only [canon] at h
    rcases min_choice (pack x) (pack (rc x)) with hx | hx <;>
    rcases min_choice (pack y) (pack (rc y)) with hy | hy <;>
    rw [hx, hy] at h
    · exact Or.inl (pack_injective hlen h).symm
    · have := pack_injective (by simp [hlen]) h
      exact Or.inr (by rw [← rc_rc y, ← this])
    · exact Or.inr (pack_injective hrc h).symm
    · have := pack_injective (by simp [hlen]) h
      exact Or.inl (by rw [← rc_rc x, ← rc_rc y, this])
  · rintro (rfl | rfl)
    · rfl
    · simp

/-! ## The word-level reverse complement -/

section Word

open BitVec

/-- `__builtin_bswap64`, written as eight byte extractions in reverse order. -/
def bswap64 (w : BitVec 64) : BitVec 64 :=
  w.extractLsb' 0 8 ++ w.extractLsb' 8 8 ++ w.extractLsb' 16 8 ++ w.extractLsb' 24 8 ++
  w.extractLsb' 32 8 ++ w.extractLsb' 40 8 ++ w.extractLsb' 48 8 ++ w.extractLsb' 56 8

/-- `reverse_complement` in kmer.hpp, step for step. -/
def rcWord (w : BitVec 64) (k : ℕ) : BitVec 64 :=
  let bits := bswap64 (~~~w)
  let bits := ((bits >>> 4) &&& 0x0F0F0F0F0F0F0F0F#64) ||| ((bits &&& 0x0F0F0F0F0F0F0F0F#64) <<< 4)
  let bits := ((bits >>> 2) &&& 0x3333333333333333#64) ||| ((bits &&& 0x3333333333333333#64) <<< 2)
  if k = 0 then 0 else bits >>> (64 - 2 * k)

/-- Base i of a packed word, counted from the low end, so the last base of the k-mer is base 0. -/
def baseAt (w : BitVec 64) (i : ℕ) : BitVec 64 := (w >>> (2 * i)) &&& 3#64

/-- The reverse complement defined base by base, where output base i is the complement of input base k - 1 - i. -/
def rcSpec (w : BitVec 64) (k : ℕ) : BitVec 64 :=
  (List.range k).foldl (fun acc i => acc ||| ((baseAt w (k - 1 - i) ^^^ 3#64) <<< (2 * i))) 0

/-- The bit trick equals the base-by-base reverse complement for every k from 0 to 32.
No mask on `w` is needed, because the final shift keeps only bits that came from the low 2k bits of the input.
`interval_cases` turns each k into a literal, after which the base loop unfolds into plain bit operations for `bv_decide`. -/
theorem rcWord_eq_rcSpec (k : ℕ) (hk : k ≤ 32) (w : BitVec 64) : rcWord w k = rcSpec w k := by
  interval_cases k <;>
    (simp only [rcWord, rcSpec, bswap64, baseAt, List.range_succ, List.range_zero, List.foldl_nil,
      List.foldl_cons, List.foldl_append, List.nil_append]; bv_decide)

end Word

end OptikProofs.Kmer
