import Std.Tactic.BVDecide

/-!
# The spill files' varint coding round-trips

The coding is `put_varint` and `save_nodes` in `src/kmer_count.cpp`, read back by the reader's `get`.
Past the memory budget, the counter writes each partition's solid k-mers to disk this way and reads them back for the graph build.
A wrong byte there would change counts or edges only on spilled runs, so the in-memory checks would never catch it.

Values are LEB128, with seven payload bits per byte, low bits first, and bit 7 set on every byte but the last.
Words are `BitVec 64` and bytes `BitVec 8`, as the code's `uint64_t` and `uint8_t`, so wrap-around is modelled too.
`decode_encode` shows that reading a value back gives the value written and leaves the bytes after it untouched.
By `decodeNodes_encodeNodes`, a whole list (key deltas, counts and one mask byte per k-mer) reads back unchanged.
The deltas need no sortedness for this, because a wrapping subtraction and a wrapping addition cancel.
The code sorts the keys only to make the deltas short.
-/

namespace OptikProofs.Varint

/-- The bytes `put_varint` writes for one value. -/
def encode (v : BitVec 64) : List (BitVec 8) :=
  if h : v < 0x80 then [v.setWidth 8]
  else ((v ||| 0x80).setWidth 8) :: encode (v >>> 7)
termination_by v.toNat
decreasing_by
  rw [BitVec.lt_def] at h
  simp only [BitVec.toNat_ushiftRight, Nat.shiftRight_eq_div_pow]
  simp at h
  omega

/-- The reader's `get`, from bit `shift` on with the bits read so far in `acc`.
It returns the value and the bytes after it, or `none` where the code throws "spill file corrupt". -/
def decodeFrom (shift : Nat) (acc : BitVec 64) : List (BitVec 8) → Option (BitVec 64 × List (BitVec 8))
  | [] => none
  | b :: rest =>
    let acc' := acc ||| ((b &&& 0x7F).setWidth 64 <<< shift)
    if b &&& 0x80 = 0 then some (acc', rest) else decodeFrom (shift + 7) acc' rest

theorem last_byte_payload (v : BitVec 64) (h : v < 0x80) : (v.setWidth 8 &&& 0x7F).setWidth 64 = v := by
  bv_decide

theorem last_byte_stops (v : BitVec 64) (h : v < 0x80) : v.setWidth 8 &&& 0x80 = 0 := by
  bv_decide

theorem more_byte_payload (v : BitVec 64) : ((v ||| 0x80).setWidth 8 &&& 0x7F).setWidth 64 = v &&& 0x7F := by
  bv_decide

theorem more_byte_continues (v : BitVec 64) : ¬((v ||| 0x80).setWidth 8 &&& 0x80 = 0) := by
  bv_decide

theorem split_low7 (v : BitVec 64) : (v &&& 0x7F) ||| ((v >>> 7) <<< 7) = v := by
  bv_decide

/-- Reading the bytes of `v` from bit `shift` ORs `v << shift` into the accumulator and stops right after them. -/
theorem decodeFrom_encode (v : BitVec 64) :
    ∀ (shift : Nat) (acc : BitVec 64) (rest : List (BitVec 8)),
      decodeFrom shift acc (encode v ++ rest) = some (acc ||| (v <<< shift), rest) := by
  induction v using encode.induct with
  | case1 v h =>
    intro shift acc rest
    rw [encode, dif_pos h]
    simp only [List.cons_append, List.nil_append, decodeFrom, last_byte_stops v h, last_byte_payload v h]
    simp
  | case2 v h ih =>
    intro shift acc rest
    rw [encode, dif_neg h]
    simp only [List.cons_append, decodeFrom, more_byte_continues v, more_byte_payload v, if_false]
    rw [ih (shift + 7)]
    have e : v <<< shift = (v &&& 0x7F) <<< shift ||| (v >>> 7) <<< (shift + 7) :=
      calc v <<< shift = ((v &&& 0x7F) ||| ((v >>> 7) <<< 7)) <<< shift := by rw [split_low7]
        _ = (v &&& 0x7F) <<< shift ||| (v >>> 7) <<< (shift + 7) := by
          rw [BitVec.shiftLeft_or_distrib, ← BitVec.shiftLeft_add, Nat.add_comm 7 shift]
          rfl
    rw [e, BitVec.or_assoc]

theorem decode_encode (v : BitVec 64) (rest : List (BitVec 8)) : decodeFrom 0 0 (encode v ++ rest) = some (v, rest) := by
  rw [decodeFrom_encode]
  simp

/-! ## Solid k-mer lists -/

/-- A spilled k-mer holds its key, its count (`uint32_t`) and its edge mask (one byte). -/
structure Node where
  key : BitVec 64
  cnt : BitVec 32
  mask : BitVec 8
  deriving DecidableEq

/-- `save_nodes` writes each k-mer's key delta from the previous key, then its count, then its mask byte. -/
def encodeNodes (prev : BitVec 64) : List Node → List (BitVec 8)
  | [] => []
  | n :: ns => encode (n.key - prev) ++ encode (n.cnt.setWidth 64) ++ [n.mask] ++ encodeNodes n.key ns

/-- The reader's loop over `ref.count` k-mers, with the running key sum. -/
def decodeNodes (key : BitVec 64) : Nat → List (BitVec 8) → Option (List Node)
  | 0, _ => some []
  | c + 1, bytes =>
    match decodeFrom 0 0 bytes with
    | none => none
    | some (delta, r1) =>
      match decodeFrom 0 0 r1 with
      | none => none
      | some (cnt, r2) =>
        match r2 with
        | [] => none
        | m :: r3 =>
          match decodeNodes (key + delta) c r3 with
          | none => none
          | some rest => some (⟨key + delta, cnt.setWidth 32, m⟩ :: rest)

theorem setWidth_roundtrip (c : BitVec 32) : (c.setWidth 64).setWidth 32 = c := by
  bv_decide

/-- A spilled list reads back exactly as it was written, whatever its order. -/
theorem decodeNodes_encodeNodes (ns : List Node) :
    ∀ prev : BitVec 64, decodeNodes prev ns.length (encodeNodes prev ns) = some ns := by
  induction ns with
  | nil => intro prev; rfl
  | cons n ns ih =>
    intro prev
    have hkey : prev + (n.key - prev) = n.key := by rw [BitVec.add_comm, BitVec.sub_add_cancel]
    rw [List.length_cons, encodeNodes, List.append_assoc, List.append_assoc, decodeNodes]
    simp only [decode_encode, List.singleton_append, hkey, ih n.key, setWidth_roundtrip]

end OptikProofs.Varint
