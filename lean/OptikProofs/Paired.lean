import Mathlib.Order.MinMax
import Mathlib.Tactic.Linarith

/-!
# Pair geometry and TLEN signs

Models `fr_tlen`, `proper` and the TLEN part of `set_pair_fields` in `tools/src/map/paired.rs`.

`fr_tlen` and `proper` do not depend on which mate is passed first, so pairing reads 1 and 2 the other way round gives the same answer.
`set_pair_fields` writes the template length on both primaries.
The SAM specification (section 1.4, field 9) requires a positive TLEN on the leftmost mate and the same value negated on the other.
`tlen_opposite` and `tlen_left_positive` show the code does exactly that.
`tlen_ne_zero_iff` shows it is non-zero exactly when both mates are mapped to one reference.
-/

namespace OptikProofs.Paired

/-- A placement's reference id, strand and half-open reference span, as `Place` holds them. -/
structure Place where
  rid : ℕ
  rev : Bool
  start : ℤ
  stop : ℤ

/-- `fr_tlen`, the FR template length from the forward mate's start to the reverse mate's end. -/
def frTlen (a b : Place) : Option ℤ :=
  if a.rid ≠ b.rid ∨ a.rev = b.rev then none
  else
    let (f, r) := if a.rev then (b, a) else (a, b)
    if f.start > r.stop then none else some (r.stop - f.start)

theorem frTlen_comm (a b : Place) : frTlen a b = frTlen b a := by
  obtain ⟨ar, av, as, ae⟩ := a
  obtain ⟨br, bv, bs, be⟩ := b
  unfold frTlen
  by_cases h : ar = br
  · subst h
    cases av <;> cases bv <;> simp
  · cases av <;> cases bv <;> simp [h, Ne.symm h]

/-- `proper`, with the insert model's bounds. -/
def proper (a b : Place) (lo hi : ℤ) : Prop := ∃ t, frTlen a b = some t ∧ lo ≤ t ∧ t ≤ hi

theorem proper_comm (a b : Place) (lo hi : ℤ) : proper a b lo hi ↔ proper b a lo hi := by
  simp [proper, frTlen_comm a b]

/-- The fields of a primary record that `set_pair_fields` reads, with `pos` and `stop` its 1-based first and last aligned base. -/
structure Prim where
  mapped : Bool
  rname : String
  pos : ℤ
  rev : Bool
  stop : ℤ

/-- The template length both primaries share, from the outer coordinates, or 0. -/
def tlen (p1 p2 : Prim) : ℤ :=
  if p1.mapped ∧ p2.mapped ∧ p1.rname = p2.rname then max p1.stop p2.stop - min p1.pos p2.pos + 1 else 0

/-- Whether mate 1 is the leftmost, with the code's tie rule. -/
def firstIsLeft (p1 p2 : Prim) : Prop := p1.pos < p2.pos ∨ (p1.pos = p2.pos ∧ (p1.rev = p2.rev ∨ p1.rev = false))

instance (p1 p2 : Prim) : Decidable (firstIsLeft p1 p2) := by unfold firstIsLeft; infer_instance

/-- The TLEN fields of mate 1's and mate 2's primaries. -/
def tlenFields (p1 p2 : Prim) : ℤ × ℤ :=
  let t := tlen p1 p2
  if firstIsLeft p1 p2 then (t, -t) else (-t, t)

theorem tlen_opposite (p1 p2 : Prim) : (tlenFields p1 p2).1 = -(tlenFields p1 p2).2 := by
  unfold tlenFields
  split_ifs <;> simp

/-- A record's span is non-empty. -/
def Prim.Wf (p : Prim) : Prop := p.pos ≤ p.stop

theorem tlen_pos {p1 p2 : Prim} (h1 : p1.Wf) (hm : p1.mapped ∧ p2.mapped ∧ p1.rname = p2.rname) :
    0 < tlen p1 p2 := by
  simp only [tlen, if_pos hm]
  have : min p1.pos p2.pos ≤ p1.pos := min_le_left _ _
  have : p1.stop ≤ max p1.stop p2.stop := le_max_left _ _
  unfold Prim.Wf at h1
  linarith

/-- The leftmost mate's TLEN is the non-negative one. -/
theorem tlen_left_positive {p1 p2 : Prim} (h1 : p1.Wf) (hleft : p1.pos < p2.pos) :
    0 ≤ (tlenFields p1 p2).1 ∧ (tlenFields p1 p2).2 ≤ 0 := by
  have hf : firstIsLeft p1 p2 := Or.inl hleft
  have h0 : 0 ≤ tlen p1 p2 := by
    by_cases hm : p1.mapped ∧ p2.mapped ∧ p1.rname = p2.rname
    · exact le_of_lt (tlen_pos h1 hm)
    · simp [tlen, if_neg hm]
  simp only [tlenFields, if_pos hf]
  constructor <;> linarith

/-- TLEN is non-zero exactly when both mates map to one reference. -/
theorem tlen_ne_zero_iff {p1 p2 : Prim} (h1 : p1.Wf) :
    tlen p1 p2 ≠ 0 ↔ p1.mapped ∧ p2.mapped ∧ p1.rname = p2.rname := by
  constructor
  · intro h
    by_contra hm
    exact h (by simp [tlen, if_neg hm])
  · intro hm
    exact ne_of_gt (tlen_pos h1 hm)

end OptikProofs.Paired
