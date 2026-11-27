import Mathlib.Data.Nat.ModEq
import Std.Tactic.BVDecide
import OptikProofs.Spss

/-!
# The minimizer order hashes are injective

`fmix64` transcribes `mix64(value, 0)` in `include/wyhash.hpp`, which `SpssCore::mhash` (`include/spss_core.hpp`) uses to order canonical m-mers.
`wang` transcribes `minimizer_order_hash(key, mask)` in `include/minimizer.hpp`, which orders the seed index's canonical k-mers.
Both are made of two kinds of step, and each kind is injective on its own.
* An xor with a right shift of itself, which `bv_decide` checks for each shift the code uses.
* A multiplication by an odd constant, possibly plus a constant, modulo 2^n.
  An odd number is coprime to 2^n, so it cancels.

`Spss.present_iff` assumes the order is injective on packed m-mers (`Spss.OrdInj`).
`mhash_ordInj` gives exactly that for the code's order, so for `mhash` the hypothesis holds whenever m ≤ 32.
-/

namespace OptikProofs.Hash

/-- One xor-shift step, `x ^= x >> s`. -/
def xs (s : ℕ) (x : BitVec 64) : BitVec 64 := x ^^^ (x >>> s)

theorem xs33_inj {x y : BitVec 64} (h : xs 33 x = xs 33 y) : x = y := by
  simp only [xs] at h
  bv_decide

theorem xs24_inj {x y : BitVec 64} (h : xs 24 x = xs 24 y) : x = y := by
  simp only [xs] at h
  bv_decide

theorem xs14_inj {x y : BitVec 64} (h : xs 14 x = xs 14 y) : x = y := by
  simp only [xs] at h
  bv_decide

theorem xs28_inj {x y : BitVec 64} (h : xs 28 x = xs 28 y) : x = y := by
  simp only [xs] at h
  bv_decide

/-- A multiplier with a 64-bit inverse cancels. -/
theorem mul_cancel {c cinv x y : BitVec 64} (hinv : cinv * c = 1) (h : c * x = c * y) : x = y := by
  have hx : x = cinv * (c * x) := by rw [← BitVec.mul_assoc, hinv]; simp
  have hy : y = cinv * (c * y) := by rw [← BitVec.mul_assoc, hinv]; simp
  rw [hx, hy, h]

/-! ## fmix64 -/

/-- MurmurHash3's two multipliers. -/
def c1 : BitVec 64 := 0xff51afd7ed558ccd
def c2 : BitVec 64 := 0xc4ceb9fe1a85ec53

/-- `mix64(value, seed)`, step for step. -/
def fmix64 (value seed : BitVec 64) : BitVec 64 :=
  xs 33 (c2 * xs 33 (c1 * xs 33 (value ^^^ seed)))

theorem c1_inv : (0x4f74430c22a54005 : BitVec 64) * c1 = 1 := by decide
theorem c2_inv : (0x9cb4b2f8129337db : BitVec 64) * c2 = 1 := by decide

/-- For any seed, `mix64` is a bijection of 64-bit words, so distinct values never collide. -/
theorem fmix64_injective (seed : BitVec 64) : Function.Injective (fun v => fmix64 v seed) := by
  intro x y h
  simp only [fmix64] at h
  have h1 := mul_cancel c1_inv (xs33_inj (mul_cancel c2_inv (xs33_inj h)))
  have h2 := xs33_inj h1
  have h3 : x ^^^ seed ^^^ seed = y ^^^ seed ^^^ seed := by rw [h2]
  simpa [BitVec.xor_assoc] using h3

/-- `mhash` on a packed m-mer, as a function of the packed value. -/
def mhash (n : ℕ) : ℕ := (fmix64 (BitVec.ofNat 64 n) 0).toNat

/-- The order `SpssCore` puts on canonical m-mers is injective on every value below 2^64, which covers m ≤ 32. -/
theorem mhash_injOn {a b : ℕ} (ha : a < 2 ^ 64) (hb : b < 2 ^ 64) (h : mhash a = mhash b) : a = b := by
  have := fmix64_injective 0 (BitVec.eq_of_toNat_eq h)
  have := congrArg BitVec.toNat this
  simp only [BitVec.toNat_ofNat] at this
  rwa [Nat.mod_eq_of_lt ha, Nat.mod_eq_of_lt hb] at this

/-- So `Spss.present_iff`'s order hypothesis holds for the code's `mhash`, at every m ≤ 32. -/
theorem mhash_ordInj {m : ℕ} (hm : m ≤ 32) : Spss.OrdInj m mhash := by
  have h4 : 4 ^ m ≤ 2 ^ 64 := by
    rw [show (4 : ℕ) = 2 ^ 2 by rfl, ← Nat.pow_mul]
    exact Nat.pow_le_pow_right (by decide) (by omega)
  intro a b ha hb h
  exact mhash_injOn (by omega) (by omega) h

/-! ## Wang's hash, masked to 2n bits -/

/-- A value fits in n bits. -/
def Fits (n : ℕ) (x : BitVec 64) : Prop := x.toNat < 2 ^ n

/-- The mask of n low bits, `kmer_mask` for n = 2k. -/
def maskOf (n : ℕ) : BitVec 64 := BitVec.ofNat 64 (2 ^ n - 1)

theorem toNat_and_mask {n : ℕ} (hn : n ≤ 64) (a : BitVec 64) : (a &&& maskOf n).toNat = a.toNat % 2 ^ n := by
  have hm : (maskOf n).toNat = 2 ^ n - 1 := by
    simp only [maskOf, BitVec.toNat_ofNat]
    apply Nat.mod_eq_of_lt
    have : 2 ^ n ≤ 2 ^ 64 := Nat.pow_le_pow_right (by decide) hn
    have : 0 < 2 ^ n := Nat.two_pow_pos n
    omega
  rw [BitVec.toNat_and, hm, Nat.and_two_pow_sub_one_eq_mod]

theorem fits_and_mask {n : ℕ} (hn : n ≤ 64) (a : BitVec 64) : Fits n (a &&& maskOf n) := by
  simp only [Fits, toNat_and_mask hn]
  exact Nat.mod_lt _ (Nat.two_pow_pos n)

theorem fits_xs {n s : ℕ} {x : BitVec 64} (hx : Fits n x) : Fits n (xs s x) := by
  simp only [Fits, xs, BitVec.toNat_xor, BitVec.toNat_ushiftRight] at hx ⊢
  exact Nat.xor_lt_two_pow hx (lt_of_le_of_lt (Nat.shiftRight_le _ _) hx)

/-- A masked step whose value is c·x + d modulo 2^64, with c odd, is injective on values of n bits. -/
theorem affine_inj {n : ℕ} (hn : n ≤ 64) (f : BitVec 64 → BitVec 64) (c d : ℕ) (hc : c % 2 = 1)
    (hf : ∀ x : BitVec 64, (f x).toNat ≡ c * x.toNat + d [MOD 2 ^ 64])
    {x y : BitVec 64} (hx : Fits n x) (hy : Fits n y) (h : f x &&& maskOf n = f y &&& maskOf n) : x = y := by
  have hdvd : 2 ^ n ∣ 2 ^ 64 := Nat.pow_dvd_pow 2 hn
  have hxy : (f x).toNat ≡ (f y).toNat [MOD 2 ^ n] := by
    have := congrArg BitVec.toNat h
    rwa [toNat_and_mask hn, toNat_and_mask hn] at this
  have hlin : c * x.toNat + d ≡ c * y.toNat + d [MOD 2 ^ n] :=
    (((hf x).of_dvd hdvd).symm.trans hxy).trans ((hf y).of_dvd hdvd)
  have h2 : Nat.Coprime 2 c := by
    unfold Nat.Coprime
    rw [Nat.gcd_rec, hc]
    rfl
  have hcop : Nat.Coprime (2 ^ n) c := Nat.Coprime.pow_left n h2
  have hmod := Nat.ModEq.cancel_left_of_coprime hcop (Nat.ModEq.add_right_cancel' d hlin)
  apply BitVec.eq_of_toNat_eq
  rw [Nat.ModEq, Nat.mod_eq_of_lt hx, Nat.mod_eq_of_lt hy] at hmod
  exact hmod

def w1 (x : BitVec 64) : BitVec 64 := ~~~x + (x <<< 21)
def w3 (x : BitVec 64) : BitVec 64 := (x + (x <<< 3)) + (x <<< 8)
def w5 (x : BitVec 64) : BitVec 64 := (x + (x <<< 2)) + (x <<< 4)
def w7 (x : BitVec 64) : BitVec 64 := x + (x <<< 31)

theorem w1_affine (x : BitVec 64) : (w1 x).toNat ≡ 2097151 * x.toNat + 18446744073709551615 [MOD 2 ^ 64] := by
  have := x.isLt
  simp only [w1, Nat.ModEq, BitVec.toNat_add, BitVec.toNat_not, BitVec.toNat_shiftLeft, Nat.shiftLeft_eq] at *
  simp only [Nat.reducePow] at *
  omega

theorem w3_affine (x : BitVec 64) : (w3 x).toNat ≡ 265 * x.toNat + 0 [MOD 2 ^ 64] := by
  have := x.isLt
  simp only [w3, Nat.ModEq, BitVec.toNat_add, BitVec.toNat_shiftLeft, Nat.shiftLeft_eq] at *
  simp only [Nat.reducePow] at *
  omega

theorem w5_affine (x : BitVec 64) : (w5 x).toNat ≡ 21 * x.toNat + 0 [MOD 2 ^ 64] := by
  have := x.isLt
  simp only [w5, Nat.ModEq, BitVec.toNat_add, BitVec.toNat_shiftLeft, Nat.shiftLeft_eq] at *
  simp only [Nat.reducePow] at *
  omega

theorem w7_affine (x : BitVec 64) : (w7 x).toNat ≡ 2147483649 * x.toNat + 0 [MOD 2 ^ 64] := by
  have := x.isLt
  simp only [w7, Nat.ModEq, BitVec.toNat_add, BitVec.toNat_shiftLeft, Nat.shiftLeft_eq] at *
  simp only [Nat.reducePow] at *
  omega

/-- `minimizer_order_hash(key, mask)`, step for step. -/
def wang (mask x : BitVec 64) : BitVec 64 :=
  let x := w1 x &&& mask
  let x := xs 24 x
  let x := w3 x &&& mask
  let x := xs 14 x
  let x := w5 x &&& mask
  let x := xs 28 x
  w7 x &&& mask

/-- With the mask of n ≤ 64 bits, the hash is injective on keys of n bits.
So two distinct k-mers (n = 2k) never tie, as the comment on `minimizer_order_hash` claims. -/
theorem wang_injOn {n : ℕ} (hn : n ≤ 64) {x y : BitVec 64} (hx : Fits n x) (hy : Fits n y)
    (h : wang (maskOf n) x = wang (maskOf n) y) : x = y := by
  simp only [wang] at h
  have f1 := fun z => fits_and_mask hn (w1 z)
  have f3 := fun z => fits_and_mask hn (w3 z)
  have f5 := fun z => fits_and_mask hn (w5 z)
  have h7 := affine_inj hn w7 _ _ (by decide) w7_affine (fits_xs (f5 _)) (fits_xs (f5 _)) h
  have h5 := affine_inj hn w5 _ _ (by decide) w5_affine (fits_xs (f3 _)) (fits_xs (f3 _)) (xs28_inj h7)
  have h3 := affine_inj hn w3 _ _ (by decide) w3_affine (fits_xs (f1 _)) (fits_xs (f1 _)) (xs14_inj h5)
  exact affine_inj hn w1 _ _ (by decide) w1_affine hx hy (xs24_inj h3)

end OptikProofs.Hash
