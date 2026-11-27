import OptikProofs.Hash
import OptikProofs.Dbg
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-!
# Rolling extraction of canonical k-mers

The loop modelled here is `extract_kmers` in `include/kmer_extractor.hpp`.
The same two updates run in the counter's loops in `src/kmer_count.cpp`.

The loop holds a forward word, a reverse-complement word and the number of bases since the last N, in 64-bit words.
Each base shifts the forward word left and puts the base at the bottom.
It shifts the reverse word right and puts the complemented base at bit 2(k - 1).
Both are then masked to 2k bits.
`run_spec` shows that every base, once k have been read since the last N, emits `canon` of the k bases ending at it.
That is the `Kmer.canon` the rest of the proofs use.

The invariant behind it says the forward word is `pack p % 4^k`, the last k bases read.
The reverse word is the first k bases of `rc p` padded with A, which also covers the first k - 1 bases, before a full k-mer is in.
-/

namespace OptikProofs.Rolling

open OptikProofs.Kmer OptikProofs.Hash

/-! ## Packing lemmas -/

theorem pack_foldl (acc : ℕ) (b : List Base) :
    b.foldl (fun acc x => acc * 4 + x.val) acc = acc * 4 ^ b.length + pack b := by
  induction b generalizing acc with
  | nil => simp [pack]
  | cons x b ih =>
    simp only [List.foldl_cons, pack, List.length_cons] at ih ⊢
    rw [ih, ih (0 * 4 + x.val)]
    ring

theorem pack_append_list (a b : List Base) : pack (a ++ b) = pack a * 4 ^ b.length + pack b := by
  rw [pack, List.foldl_append]
  exact pack_foldl _ b

theorem pack_cons (c : Base) (t : List Base) : pack (c :: t) = c.val * 4 ^ t.length + pack t := by
  have := pack_append_list [c] t
  simpa [pack] using this

theorem pack_replicate_zero (n : ℕ) : pack (List.replicate n 0) = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [List.replicate_succ, pack_cons, ih]; simp

/-- Masking a packed word to its low 2j bits leaves its last j bases. -/
theorem pack_mod_eq (l : List Base) (j : ℕ) : pack l % 4 ^ j = pack (l.drop (l.length - j)) := by
  by_cases hj : l.length ≤ j
  · rw [show l.length - j = 0 by omega, List.drop_zero]
    exact Nat.mod_eq_of_lt (lt_of_lt_of_le (pack_lt l) (Nat.pow_le_pow_right (by decide) hj))
  · have hsplit := List.take_append_drop (l.length - j) l
    have hd : (l.drop (l.length - j)).length = j := by simp; omega
    conv_lhs => rw [← hsplit, pack_append_list, hd]
    have hlt := pack_lt (l.drop (l.length - j))
    rw [hd] at hlt
    rw [Nat.add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt hlt]

theorem rc_append_single (p : List Base) (b : Base) : rc (p ++ [b]) = comp b :: rc p := by
  simp [rc]

/-! ## The reverse word's value -/

/-- The first k bases of `rc p`, padded with A to k bases. -/
def revSpec (k : ℕ) (p : List Base) : ℕ := pack ((rc p ++ List.replicate k 0).take k)

theorem revSpec_lt (k : ℕ) (p : List Base) : revSpec k p < 4 ^ k := by
  have h := pack_lt ((rc p ++ List.replicate k 0).take k)
  have hl : ((rc p ++ List.replicate k 0).take k).length = k := by simp
  rwa [hl] at h

theorem revSpec_step {k : ℕ} (hk : 1 ≤ k) (p : List Base) (b : Base) :
    revSpec k (p ++ [b]) = revSpec k p / 4 + (comp b).val * 4 ^ (k - 1) := by
  set L := rc p ++ List.replicate k 0 with hL
  have hlen : k ≤ L.length := by simp [hL]
  have htake : L.take k = L.take (k - 1) ++ [L[k - 1]'(by omega)] := by
    conv_lhs => rw [show k = k - 1 + 1 by omega]
    rw [List.take_add_one, List.getElem?_eq_getElem (by omega)]
    rfl
  have hlow : pack (L.take (k - 1)) = pack (L.take k) / 4 := by
    rw [htake, pack_append]
    have := (L[k - 1]'(by omega)).isLt
    omega
  have hl1 : (L.take (k - 1)).length = k - 1 := by simp; omega
  unfold revSpec
  rw [rc_append_single, List.cons_append, show k = k - 1 + 1 by omega, List.take_succ_cons, pack_cons,
    show k - 1 + 1 = k by omega, ← hL, hl1, hlow]
  ring

/-- Once k bases are in, the reverse word is the packed reverse complement of the last k. -/
theorem revSpec_full {k : ℕ} {p : List Base} (hk : k ≤ p.length) :
    revSpec k p = pack (rc (p.drop (p.length - k))) := by
  rw [revSpec, List.take_append_of_le_length (by simpa using hk), Dbg.take_rc]

/-! ## The code's word updates -/

theorem two_pow_two_mul (k : ℕ) : 2 ^ (2 * k) = 4 ^ k := by
  rw [Nat.pow_mul]

theorem or_low_eq_add (x bb : BitVec 64) (hb : bb < 4) : (x <<< 2) ||| bb = x * 4 + bb := by
  bv_decide

/-- `fwd = ((fwd << 2) | base) & mask` keeps the last k bases. -/
theorem stepF_toNat {k : ℕ} (hk : k ≤ 32) (x : BitVec 64) (b : Base) :
    (((x <<< 2) ||| BitVec.ofNat 64 b.val) &&& maskOf (2 * k)).toNat = (x.toNat * 4 + b.val) % 4 ^ k := by
  have hb := b.isLt
  have hbv : (BitVec.ofNat 64 b.val).toNat = b.val := by
    rw [BitVec.toNat_ofNat]
    exact Nat.mod_eq_of_lt (by omega)
  have hlt : BitVec.ofNat 64 b.val < 4 := by
    rw [BitVec.lt_def, hbv]
    exact hb
  have h4 : (4 : BitVec 64).toNat = 4 := rfl
  rw [or_low_eq_add _ _ hlt, toNat_and_mask (by omega), BitVec.toNat_add, BitVec.toNat_mul, h4, hbv,
    two_pow_two_mul]
  have hdvd : 4 ^ k ∣ 2 ^ 64 := by
    rw [← two_pow_two_mul]
    exact Nat.pow_dvd_pow 2 (by omega)
  rw [Nat.mod_mod_of_dvd _ hdvd]
  exact (((Nat.mod_modEq _ _).of_dvd hdvd).add_right _)

theorem xor3_eq_comp (b : Base) : BitVec.ofNat 64 b.val ^^^ 3 = BitVec.ofNat 64 (comp b).val := by
  revert b
  decide

/-- `rev = ((rev >> 2) | (comp << 2(k-1))) & mask` drops the last base of the reverse word and puts the new one first. -/
theorem stepR_toNat {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 32) (r : BitVec 64) (hr : r.toNat < 4 ^ k) (b : Base) :
    (((r >>> 2) ||| ((BitVec.ofNat 64 b.val ^^^ 3) <<< (2 * (k - 1)))) &&& maskOf (2 * k)).toNat =
      r.toNat / 4 + (comp b).val * 4 ^ (k - 1) := by
  set s := 2 * (k - 1) with hs
  set c := (comp b).val with hc
  have hc4 : c < 4 := (comp b).isLt
  have hP : 4 * 2 ^ s = 4 ^ k := by
    rw [hs, two_pow_two_mul, ← Nat.pow_succ']
    congr 1
    omega
  have hPle : 4 * 2 ^ s ≤ 2 ^ 64 := by
    rw [show (4 : ℕ) = 2 ^ 2 by rfl, ← Nat.pow_add]
    exact Nat.pow_le_pow_right (by decide) (by omega)
  have hrlow : r.toNat / 2 ^ 2 < 2 ^ s := Nat.div_lt_of_lt_mul (by rw [show (2 : ℕ) ^ 2 = 4 from rfl, hP]; exact hr)
  have hcs : c * 2 ^ s ≤ 3 * 2 ^ s := Nat.mul_le_mul_right _ (by omega)
  have hshift : ((BitVec.ofNat 64 c) <<< s).toNat = c * 2 ^ s := by
    have hc' : c % 2 ^ 64 = c := Nat.mod_eq_of_lt (by omega)
    rw [BitVec.toNat_shiftLeft, BitVec.toNat_ofNat, hc', Nat.shiftLeft_eq]
    exact Nat.mod_eq_of_lt (by omega)
  have hor : (r >>> 2 ||| (BitVec.ofNat 64 c) <<< s).toNat = r.toNat / 4 + c * 2 ^ s := by
    rw [BitVec.toNat_or, hshift, BitVec.toNat_ushiftRight, Nat.shiftRight_eq_div_pow, Nat.or_comm,
      ← Nat.shiftLeft_eq, ← Nat.shiftLeft_add_eq_or_of_lt hrlow, Nat.shiftLeft_eq]
    norm_num
    omega
  have hlow4 : r.toNat / 4 < 2 ^ s := by simpa using hrlow
  have hP' : 2 ^ s = 4 ^ (k - 1) := by rw [hs, two_pow_two_mul]
  rw [hP'] at hor hcs hlow4 hP
  rw [xor3_eq_comp, ← hc, toNat_and_mask (by omega), hor, two_pow_two_mul, Nat.mod_eq_of_lt (by omega)]

/-! ## The loop -/

/-- The loop state holds `fwd`, `rev` and `valid`, the count of bases since the last restart. -/
structure St where
  fwd : BitVec 64
  rev : BitVec 64
  valid : ℕ

def init : St := ⟨0, 0, 0⟩

/-- The loop body's update for one ACGT base. -/
def step (k : ℕ) (st : St) (b : Base) : St :=
  ⟨((st.fwd <<< 2) ||| BitVec.ofNat 64 b.val) &&& maskOf (2 * k),
   ((st.rev >>> 2) ||| ((BitVec.ofNat 64 b.val ^^^ 3) <<< (2 * (k - 1)))) &&& maskOf (2 * k),
   st.valid + 1⟩

/-- `can = (fwd < rev) ? fwd : rev`, passed to the callback once `valid >= kmer_len`. -/
def emit (k : ℕ) (st : St) : List ℕ :=
  if k ≤ st.valid then [(if st.fwd < st.rev then st.fwd else st.rev).toNat] else []

/-- The whole loop, with `none` for an N or any other non-base, which restarts both words. -/
def run (k : ℕ) : St → List (Option Base) → List ℕ
  | _, [] => []
  | _, none :: rest => run k init rest
  | st, some b :: rest => emit k (step k st b) ++ run k (step k st b) rest

/-- What the loop should emit after the bases `p`, for the bases `s` still to come. -/
def spec (k : ℕ) (p : List Base) : List Base → List ℕ
  | [] => []
  | b :: s => (if k ≤ p.length + 1 then [canon ((p ++ [b]).drop (p.length + 1 - k))] else []) ++ spec k (p ++ [b]) s

/-- The state after reading the bases `p` since the last restart. -/
def Inv (k : ℕ) (st : St) (p : List Base) : Prop :=
  st.valid = p.length ∧ st.fwd.toNat = pack p % 4 ^ k ∧ st.rev.toNat = revSpec k p

theorem inv_init (k : ℕ) : Inv k init [] := by
  refine ⟨rfl, by simp [init, pack], ?_⟩
  simp [init, revSpec, rc, pack_replicate_zero]

theorem inv_step {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 32) {st : St} {p : List Base} (h : Inv k st p) (b : Base) :
    Inv k (step k st b) (p ++ [b]) := by
  obtain ⟨hv, hf, hr⟩ := h
  refine ⟨by simp [step, hv], ?_, ?_⟩
  · simp only [step]
    rw [stepF_toNat hk, hf, pack_append]
    exact (((Nat.mod_modEq _ _).mul_right 4).add_right _)
  · simp only [step]
    rw [stepR_toNat hk1 hk _ (hr ▸ revSpec_lt k p), hr, revSpec_step hk1]

theorem emit_val {k : ℕ} {st : St} {p : List Base} (h : Inv k st p) (hk : k ≤ p.length) :
    (if st.fwd < st.rev then st.fwd else st.rev).toNat = canon (p.drop (p.length - k)) := by
  obtain ⟨_, hf, hr⟩ := h
  rw [pack_mod_eq] at hf
  rw [revSpec_full hk] at hr
  simp only [canon]
  split_ifs with hlt
  · rw [BitVec.lt_def] at hlt
    omega
  · rw [BitVec.lt_def] at hlt
    omega

/-- The loop emits what `spec` says, from any state that matches the bases read so far. -/
theorem run_from {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 32) (s : List Base) :
    ∀ (st : St) (p : List Base), Inv k st p → run k st (s.map some) = spec k p s := by
  induction s with
  | nil => intro st p _; rfl
  | cons b s ih =>
    intro st p h
    have h' := inv_step hk1 hk h b
    simp only [List.map_cons, run, spec]
    rw [ih _ _ h']
    congr 1
    unfold emit
    rw [h'.1]
    simp only [List.length_append, List.length_singleton]
    by_cases hle : k ≤ p.length + 1
    · rw [if_pos hle, if_pos hle, emit_val h' (by simpa using hle)]
      simp
    · rw [if_neg hle, if_neg hle]

/-- From a restart, an N-free run of bases gives `canon` of each k-mer in order, for 1 ≤ k ≤ 32. -/
theorem run_spec {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 32) (s : List Base) : run k init (s.map some) = spec k [] s :=
  run_from hk1 hk s init [] (inv_init k)

/-- An N restarts the loop, so the output for a whole sequence is the outputs of its N-free runs, one after another. -/
theorem run_none (k : ℕ) (st : St) (rest : List (Option Base)) : run k st (none :: rest) = run k init rest := rfl

end OptikProofs.Rolling
