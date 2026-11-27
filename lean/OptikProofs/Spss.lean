import OptikProofs.Minimizer

/-!
# SPSS index exactness (prop:spss)

The code is the lookup in `include/spss_core.hpp` and `include/spss_scan.hpp`.

The index stores sequences, and for every k-mer it files the position q of the k-mer's minimizer occurrence under that minimizer's canonical form.
A query x tries each offset p where its minimum order value occurs.
For each record q in the bucket it checks one candidate window, starting at q - p when x is stored forward and at q - (k - m - p) when it is stored reverse complemented.
`present_iff` states that x is reported present exactly when x or its reverse complement occurs in some sequence.

The model makes three simplifications, each stated as an assumption below.
* A record is a (sequence, offset) pair, and every candidate is bounds-checked.
  The code checks bounds only for records flagged near a sequence end, which matters for speed but not for the answer.
* Heavy buckets look the k-mer up in a skew table instead of scanning the bucket.
  That is modelled as any candidate function returning every bucket record that verifies the query, which is the one property the proof uses.
* The MPHF is assumed only to be injective on the minimizers the index was built on, so a query whose minimizer is new may land in an unrelated bucket.
-/

namespace OptikProofs.Spss

open OptikProofs.Kmer OptikProofs.Minimizer

variable {α : Type*} [LinearOrder α] (k m : ℕ) (ord : ℕ → α)

/-- The m-mer at offset o of a k-mer. -/
def mmerAt (y : List Base) (o : ℕ) : List Base := (y.drop o).take m

/-- The order is injective on packed m-mers, the only values it is ever applied to.
`Hash.mhash_ordInj` proves this for the code's order, MurmurHash3 fmix64, when m ≤ 32. -/
def OrdInj : Prop := ∀ a b, a < 4 ^ m → b < 4 ^ m → ord a = ord b → a = b

theorem canon_mmerAt_lt (y : List Base) (o : ℕ) : canon (mmerAt m y o) < 4 ^ m :=
  lt_of_le_of_lt (min_le_left _ _)
    (lt_of_lt_of_le (pack_lt _) (Nat.pow_le_pow_right (by decide) (by simp [mmerAt])))

/-- The order value of the m-mer at offset o, a hash of its canonical form in the code. -/
def ordAt (y : List Base) (o : ℕ) : α := ord (canon (mmerAt m y o))

/-- The leftmost offset with the smallest order value, which the build uses for each stored k-mer. -/
def leftTie (y : List Base) : ℕ := scanMin (ordAt m ord y) 0 (k - m)

/-- The bucket key of a k-mer, the canonical form of its minimizer. -/
def key (y : List Base) : ℕ := canon (mmerAt m y (leftTie k m ord y))

/-- The k-mer starting at offset s of a sequence. -/
def window (str : List Base) (s : ℕ) : List Base := (str.drop s).take k

/-- x occurs in the indexed sequences. -/
def Occurs (strs : List (List Base)) (x : List Base) : Prop :=
  ∃ sid, ∃ hs : sid < strs.length, ∃ s, s + k ≤ (strs[sid]).length ∧ window k (strs[sid]) s = x

/-- The records the build files under bucket `K`, one per stored k-mer with that key, giving its sequence and the offset of its minimizer occurrence. -/
def InBucket (strs : List (List Base)) (K : ℕ) (r : ℕ × ℕ) : Prop :=
  ∃ hs : r.1 < strs.length, ∃ i, i + k ≤ (strs[r.1]).length ∧
    key k m ord (window k (strs[r.1]) i) = K ∧ r.2 = i + leftTie k m ord (window k (strs[r.1]) i)

/-- Record r confirms x at tied offset p, read forward or reverse complemented. -/
def Verifies (strs : List (List Base)) (x : List Base) (p : ℕ) (r : ℕ × ℕ) : Prop :=
  ∃ hs : r.1 < strs.length,
    (p ≤ r.2 ∧ r.2 - p + k ≤ (strs[r.1]).length ∧ window k (strs[r.1]) (r.2 - p) = x) ∨
    (k - m - p ≤ r.2 ∧ r.2 - (k - m - p) + k ≤ (strs[r.1]).length ∧
      window k (strs[r.1]) (r.2 - (k - m - p)) = rc x)

/-- Lookup with candidate records `cand x` succeeds when some tied offset of x and some candidate record confirm it. -/
def Present (strs : List (List Base)) (cand : List Base → List (ℕ × ℕ)) (x : List Base) : Prop :=
  ∃ p, Tied (ordAt m ord x) (k - m) p ∧ ∃ r ∈ cand x, Verifies k m strs x p r

variable {k m ord}

theorem window_length {str : List Base} {s : ℕ} (h : s + k ≤ str.length) : (window k str s).length = k := by
  simp [window]; omega

theorem leftTie_tied (y : List Base) : Tied (ordAt m ord y) (k - m) (leftTie k m ord y) := by
  obtain ⟨-, hle, hmin, -⟩ := scanMin_spec (ordAt m ord y) (Nat.zero_le (k - m))
  exact ⟨hle, fun o ho => hmin o (Nat.zero_le _) ho⟩

omit [LinearOrder α] in
/-- The order values of rc y, read back to front, are those of y. -/
theorem ordAt_rc {y : List Base} (hy : y.length = k) (hmk : m ≤ k) {o : ℕ} (ho : o ≤ k - m) :
    ordAt m ord (rc y) o = ordAt m ord y (k - m - o) := by
  simp only [ordAt, mmerAt]
  have := order_mirror ord y (m := m) (o := o) (by omega) (by omega)
  rw [hy] at this
  exact this

/-- Offsets of rc y are tied exactly when their mirror images are tied in y. -/
theorem tied_rc {y : List Base} (hy : y.length = k) (hmk : m ≤ k) {o : ℕ} (ho : o ≤ k - m) :
    Tied (ordAt m ord (rc y)) (k - m) o ↔ Tied (ordAt m ord y) (k - m) (k - m - o) := by
  rw [← tied_mirror (ordAt m ord y) ho]
  constructor
  · rintro ⟨ho', hmin⟩
    refine ⟨ho', fun o' ho'' => ?_⟩
    show ordAt m ord y (k - m - o) ≤ ordAt m ord y (k - m - o')
    rw [← ordAt_rc hy hmk ho, ← ordAt_rc hy hmk ho'']
    exact hmin o' ho''
  · rintro ⟨ho', hmin⟩
    refine ⟨ho', fun o' ho'' => ?_⟩
    rw [ordAt_rc hy hmk ho, ordAt_rc hy hmk ho'']
    exact hmin o' ho''

/-- With the order injective on canonical m-mers, all tied offsets of a k-mer hold the same canonical m-mer. -/
theorem canon_of_tied (hord : OrdInj m ord) {y : List Base} {o o' : ℕ}
    (ho : Tied (ordAt m ord y) (k - m) o) (ho' : Tied (ordAt m ord y) (k - m) o') :
    canon (mmerAt m y o) = canon (mmerAt m y o') :=
  hord _ _ (canon_mmerAt_lt m y o) (canon_mmerAt_lt m y o') (le_antisymm (ho.2 o' ho'.1) (ho'.2 o ho.1))

/-- A k-mer and its reverse complement have the same bucket key, so both strands look in one bucket. -/
theorem key_rc (hord : OrdInj m ord) (hmk : m ≤ k) {y : List Base} (hy : y.length = k) :
    key k m ord (rc y) = key k m ord y := by
  have ht := leftTie_tied (k := k) (m := m) (ord := ord) (rc y)
  have hle := ht.1
  have hmirror := (tied_rc hy hmk hle).mp ht
  have hmm : ((rc y).drop (leftTie k m ord (rc y))).take m =
      rc ((y.drop (y.length - m - leftTie k m ord (rc y))).take m) := mmer_rc (by rw [hy]; omega)
  simp only [key, mmerAt]
  rw [hmm, canon_rc, hy]
  exact canon_of_tied hord hmirror (leftTie_tied (k := k) (m := m) (ord := ord) y)

variable (k m ord)

/-- prop:spss.
With every bucket record that verifies x among its candidates, lookup reports x present exactly when x or rc x occurs.
`hcand` is the only property the proof needs of light buckets (the whole bucket) and of heavy buckets (the skew table's record for x). -/
theorem present_iff (hord : OrdInj m ord) (hmk : m ≤ k) (strs : List (List Base))
    (cand : List Base → List (ℕ × ℕ))
    (hcand : ∀ x r p, InBucket k m ord strs (key k m ord x) r → Verifies k m strs x p r → r ∈ cand x)
    {x : List Base} (hx : x.length = k) :
    Present k m ord strs cand x ↔ Occurs k strs x ∨ Occurs k strs (rc x) := by
  constructor
  · -- For soundness, a confirmed candidate is a real window holding x or rc x.
    rintro ⟨p, -, r, -, hs, ⟨-, hb, hw⟩ | ⟨-, hb, hw⟩⟩
    · exact Or.inl ⟨r.1, hs, _, hb, hw⟩
    · exact Or.inr ⟨r.1, hs, _, hb, hw⟩
  · rintro (⟨sid, hs, i, hb, hw⟩ | ⟨sid, hs, i, hb, hw⟩)
    · -- In the forward case the build filed (sid, i + leftTie x), and the leftmost tie is one of the offsets tried.
      set p := leftTie k m ord x
      have hp := leftTie_tied (k := k) (m := m) (ord := ord) x
      have hrec : InBucket k m ord strs (key k m ord x) (sid, i + p) := ⟨hs, i, hb, by rw [hw], by rw [hw]⟩
      have hver : Verifies k m strs x p (sid, i + p) :=
        ⟨hs, Or.inl ⟨by simp, by simpa using hb, by simpa using hw⟩⟩
      exact ⟨p, hp, _, hcand x _ p hrec hver, hver⟩
    · -- In the reverse case rc x was filed with its leftmost tie p', whose mirror k - m - p' is tied in x.
      set y := rc x
      have hy : y.length = k := by simp [y, hx]
      set p' := leftTie k m ord y
      have hp' := leftTie_tied (k := k) (m := m) (ord := ord) y
      have hle : p' ≤ k - m := hp'.1
      have hxy : rc y = x := rc_rc x
      have htied : Tied (ordAt m ord x) (k - m) (k - m - p') := by
        have := (tied_rc (ord := ord) hy hmk (show k - m - p' ≤ k - m by omega)).mpr
        rw [hxy, show k - m - (k - m - p') = p' by omega] at this
        exact this hp'
      have hkey : key k m ord y = key k m ord x := by
        rw [← key_rc hord hmk hy, hxy]
      have hrec : InBucket k m ord strs (key k m ord x) (sid, i + p') :=
        ⟨hs, i, hb, by rw [hw, hkey], by rw [hw]⟩
      have hver : Verifies k m strs x (k - m - p') (sid, i + p') := by
        refine ⟨hs, Or.inr ⟨by simp; omega, ?_, ?_⟩⟩
        · simpa [show k - m - (k - m - p') = p' by omega] using hb
        · simpa [show k - m - (k - m - p') = p' by omega] using hw
      exact ⟨k - m - p', htied, _, hcand x _ _ hrec hver, hver⟩

/-- In streaming extension, when the next query k-mer matches the stored window one step on, it is reported present without a lookup.
That shortcut never changes an answer, because such a match already makes the lookup succeed. -/
theorem extension_sound (hord : OrdInj m ord) (hmk : m ≤ k) (strs : List (List Base))
    (cand : List Base → List (ℕ × ℕ))
    (hcand : ∀ x r p, InBucket k m ord strs (key k m ord x) r → Verifies k m strs x p r → r ∈ cand x)
    {x : List Base} (hx : x.length = k) (hext : Occurs k strs x ∨ Occurs k strs (rc x)) :
    Present k m ord strs cand x :=
  (present_iff k m ord hord hmk strs cand hcand hx).mpr hext

/-- `Occurs` as a computation, which the test vectors use for their expected answers. -/
def occursB (strs : List (List Base)) (x : List Base) : Bool :=
  strs.any fun str => (List.range (str.length + 1 - k)).any fun s => window k str s == x

theorem occursB_iff (strs : List (List Base)) (x : List Base) : occursB k strs x = true ↔ Occurs k strs x := by
  simp only [occursB, List.any_eq_true, List.mem_range, beq_iff_eq, Occurs]
  constructor
  · rintro ⟨str, hstr, s, hs, hw⟩
    obtain ⟨sid, hsid, rfl⟩ := List.getElem_of_mem hstr
    exact ⟨sid, hsid, s, by omega, hw⟩
  · rintro ⟨sid, hsid, s, hs, hw⟩
    exact ⟨_, List.getElem_mem hsid, s, by omega, hw⟩

/-! ## Dense k-mer ids -/

section DenseIds

/-- k-mers held by strings of lengths `lens`, each string holding len - k + 1 of them. -/
def total (lens : List ℕ) : ℕ := (lens.map fun len => len - k + 1).sum

/-- Every k-mer position as (string index, offset in that string), in storage order, numbering strings from `i`. -/
def positionsFrom : ℕ → List ℕ → List (ℕ × ℕ)
  | _, [] => []
  | i, len :: rest => (List.range (len - k + 1)).map (fun o => (i, o)) ++ positionsFrom (i + 1) rest

/-- In `kmer_id` the k-mer at offset o of string i sits at global offset start_i + o and gets s - (k - 1) · i. -/
def idOf (lens : List ℕ) (pos : ℕ × ℕ) : ℕ := (lens.take pos.1).sum + pos.2 - (k - 1) * pos.1

theorem sum_eq_total (hk : 1 ≤ k) {lens : List ℕ} (hlens : ∀ len ∈ lens, k ≤ len) :
    lens.sum = total k lens + (k - 1) * lens.length := by
  induction lens with
  | nil => simp [total]
  | cons len rest ih =>
    have hlen := hlens len (by simp)
    have := ih fun l hl => hlens l (by simp [hl])
    simp only [total, List.map_cons, List.sum_cons, List.length_cons] at this ⊢
    rw [this, Nat.mul_succ]
    omega

theorem ids_from (hk : 1 ≤ k) (pre : List ℕ) : ∀ (rest : List ℕ), (∀ len ∈ pre ++ rest, k ≤ len) →
    (positionsFrom k pre.length rest).map (idOf k (pre ++ rest)) =
      List.range' (total k pre) (total k rest)
  | [], _ => by simp [positionsFrom, total]
  | len :: rest, hlens => by
    have hpre := sum_eq_total k hk (lens := pre) fun l hl => hlens l (by simp [hl])
    have hlen : k ≤ len := hlens len (by simp)
    have ih := ids_from hk (pre ++ [len]) rest (by simpa using hlens)
    simp only [List.length_append, List.length_singleton, List.append_assoc, List.singleton_append] at ih
    simp only [positionsFrom, List.map_append, List.map_map, Function.comp_def, ih]
    have hfirst : (List.range (len - k + 1)).map (fun o => idOf k (pre ++ len :: rest) (pre.length, o)) =
        List.range' (total k pre) (len - k + 1) := by
      rw [List.range'_eq_map_range]
      refine List.map_congr_left fun o _ => ?_
      simp only [idOf, List.take_left' rfl]
      omega
    rw [hfirst, show total k (pre ++ [len]) = total k pre + (len - k + 1) by simp [total],
      List.range'_append_1, show total k (len :: rest) = len - k + 1 + total k rest by simp [total]]

/-- Dense ids, listed in storage order, are exactly 0, 1, …, n - 1, so `kmer_id` is a bijection onto [0, n).
Every string is at least k long, because SpssCore drops shorter runs. -/
theorem ids_dense (hk : 1 ≤ k) (lens : List ℕ) (hlens : ∀ len ∈ lens, k ≤ len) :
    (positionsFrom k 0 lens).map (idOf k lens) = List.range (total k lens) := by
  have := ids_from k hk [] lens (by simpa using hlens)
  simpa [total, List.range_eq_range'] using this

end DenseIds

end OptikProofs.Spss
