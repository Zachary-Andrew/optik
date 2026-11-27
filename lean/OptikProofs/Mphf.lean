import Mathlib.Data.Nat.Count
import Mathlib.Data.Finset.Card
import Mathlib.Order.Interval.Finset.Nat

/-!
# The MPHF's slot remapping is a bijection

Models the end of `PilotMphf::try_build` and the lookup `PilotMphf::operator()` in `include/phobic_mphf.hpp`, and the offsets of `Mphf`.

The pilot search places the n keys on slots of a table of size T ≥ n.
It accepts a pilot only when every key of the bucket lands on a slot that is free and not shared with another key of that bucket.
`fits_nodup` proves that this makes the placement injective.
The keys that landed at or past n are then sent to the free slots below n, called holes, in increasing order.
By `lookup_bijective` the result is a bijection from the keys onto [0, n), which is the property `Spss.present_iff` assumes of the MPHF.
The remap table is non-decreasing (`remap_sorted`), which `EfSeq` needs, and every entry is below n (`remap_lt`), which bounds its high bits.
The same bound means the hole search never reads past slot n.
Finally, `global_injective` states that `Mphf`'s per-partition offsets separate the values of different partitions.

Slots are natural numbers and the `taken` bit vector is the set `S` of slots in use.
-/

namespace OptikProofs.Mphf

open Finset

/-! ## Placing buckets -/

/-- The pilot search's acceptance test for each bucket in turn, given the slots taken before it. -/
def Fits : List ℕ → List (List ℕ) → Prop
  | _, [] => True
  | taken, b :: bs => b.Nodup ∧ (∀ s ∈ b, s ∉ taken) ∧ Fits (taken ++ b) bs

/-- A build that placed every bucket never put two keys on one slot. -/
theorem fits_nodup : ∀ (bs : List (List ℕ)) (taken : List ℕ), taken.Nodup → Fits taken bs → (taken ++ bs.flatten).Nodup
  | [], taken, h, _ => by simpa using h
  | b :: bs, taken, h, ⟨hb, hfree, hrest⟩ => by
    have h' : (taken ++ b).Nodup := List.nodup_append.mpr ⟨h, hb, fun x hx y hy hxy => hfree y hy (hxy ▸ hx)⟩
    simpa [List.append_assoc] using fits_nodup bs (taken ++ b) h' hrest

/-! ## The remap loop -/

variable (S : Finset ℕ)

theorem exists_free (h : ℕ) : ∃ x, h ≤ x ∧ x ∉ S :=
  ⟨h + S.sup id + 1, by omega, fun hx => by have := Finset.le_sup (f := id) hx; simp only [id_eq] at this; omega⟩

/-- The first free slot at or after `h`, as `while (is_taken(hole)) ++hole` finds it. -/
noncomputable def firstFree (h : ℕ) : ℕ := Nat.find (exists_free S h)

theorem firstFree_spec (h : ℕ) : h ≤ firstFree S h ∧ firstFree S h ∉ S := Nat.find_spec (exists_free S h)

/-- The free slots below x. -/
def free (x : ℕ) : ℕ := Nat.count (· ∉ S) x

/-- The slots skipped by the search are all taken, so no free slot lies between `h` and the one found. -/
theorem free_firstFree (h : ℕ) : free S (firstFree S h) = free S h := by
  have key : ∀ x, h ≤ x → x ≤ firstFree S h → free S x = free S h := by
    intro x hx hle
    induction x, hx using Nat.le_induction with
    | base => rfl
    | succ x hx ih =>
      have hlt : x < firstFree S h := by omega
      have hin : x ∈ S := by
        by_contra hn
        exact Nat.find_min (exists_free S h) hlt ⟨hx, hn⟩
      unfold free at ih ⊢
      rw [Nat.count_succ, ih (by omega)]
      simp [hin]
  exact key _ (firstFree_spec S h).1 le_rfl

theorem free_succ_of_free {x : ℕ} (hx : x ∉ S) : free S (x + 1) = free S x + 1 := by
  unfold free; rw [Nat.count_succ]; simp [hx]

/-- The loop over slots q = n .. T - 1, from `hole` and `last`, giving one remap entry per slot. -/
noncomputable def run : List ℕ → ℕ → ℕ → List ℕ
  | [], _, _ => []
  | q :: qs, hole, last =>
    if q ∈ S then firstFree S hole :: run qs (firstFree S hole + 1) (firstFree S hole)
    else last :: run qs hole last

/-- `remap_values`, which `EfSeq` stores. -/
noncomputable def remap (n T : ℕ) : List ℕ := run S (List.range' n (T - n)) 0 0

/-- The lookup returns a slot below n as it is and sends any other slot through the remap table. -/
noncomputable def lookup (n T q : ℕ) : ℕ := if q < n then q else (remap S n T).getD (q - n) 0

theorem run_ge : ∀ (qs : List ℕ) (hole last : ℕ), last ≤ hole → ∀ y ∈ run S qs hole last, last ≤ y
  | [], _, _, _ => by simp [run]
  | q :: qs, hole, last, hl => by
    have hf := (firstFree_spec S hole).1
    unfold run
    split_ifs
    · intro y hy
      rcases List.mem_cons.mp hy with rfl | hy
      · omega
      · have := run_ge qs _ _ (Nat.le_succ _) y hy; omega
    · intro y hy
      rcases List.mem_cons.mp hy with rfl | hy
      · exact le_rfl
      · exact run_ge qs _ _ hl y hy

/-- The remap table is non-decreasing. -/
theorem run_sorted : ∀ (qs : List ℕ) (hole last : ℕ), last ≤ hole → (run S qs hole last).Pairwise (· ≤ ·)
  | [], _, _, _ => by simp [run]
  | q :: qs, hole, last, hl => by
    unfold run
    split_ifs
    · exact List.pairwise_cons.mpr ⟨run_ge S qs _ _ (Nat.le_succ _), run_sorted qs _ _ (Nat.le_succ _)⟩
    · exact List.pairwise_cons.mpr ⟨run_ge S qs _ _ hl, run_sorted qs _ _ hl⟩

theorem remap_sorted (n T : ℕ) : (remap S n T).Pairwise (· ≤ ·) := run_sorted S _ 0 0 le_rfl

/-- Every entry is below `B` when the free slots below `B` outnumber the holes still to be handed out. -/
theorem run_lt (B : ℕ) : ∀ (qs : List ℕ) (hole last : ℕ), last < B →
    free S hole + qs.countP (· ∈ S) ≤ free S B → ∀ y ∈ run S qs hole last, y < B
  | [], _, _, _, _ => by simp [run]
  | q :: qs, hole, last, hl, hc => by
    obtain ⟨-, hfree⟩ := firstFree_spec S hole
    have hs := free_succ_of_free S hfree
    have hff := free_firstFree S hole
    unfold free at hs hff hc
    unfold run
    split_ifs with hq
    · rw [List.countP_cons_of_pos (by simpa using hq)] at hc
      have hlt : firstFree S hole < B := by
        by_contra hge
        have := Nat.count_monotone (· ∉ S) (show B ≤ firstFree S hole by omega)
        omega
      intro y hy
      rcases List.mem_cons.mp hy with rfl | hy
      · exact hlt
      · exact run_lt B qs _ _ hlt (by unfold free; omega) y hy
    · rw [List.countP_cons_of_neg (by simpa using hq)] at hc
      intro y hy
      rcases List.mem_cons.mp hy with rfl | hy
      · exact hl
      · exact run_lt B qs _ _ hl (by unfold free; omega) y hy

/-- A taken slot is sent to a free slot at or after the current hole. -/
theorem run_taken : ∀ (qs : List ℕ) (hole last i q : ℕ), qs[i]? = some q → q ∈ S →
    ∃ x, (run S qs hole last)[i]? = some x ∧ x ∉ S ∧ hole ≤ x
  | [], _, _, _, _, h, _ => by simp at h
  | q' :: qs, hole, last, 0, q, h, hq => by
    simp only [List.getElem?_cons_zero, Option.some.injEq] at h
    subst h
    simp only [run, if_pos hq, List.getElem?_cons_zero, Option.some.injEq, exists_eq_left']
    exact ⟨(firstFree_spec S hole).2, (firstFree_spec S hole).1⟩
  | q' :: qs, hole, last, i + 1, q, h, hq => by
    simp only [List.getElem?_cons_succ] at h
    unfold run
    split_ifs
    · obtain ⟨x, hx, hxs, hle⟩ := run_taken qs _ _ i q h hq
      exact ⟨x, by simpa using hx, hxs, by have := (firstFree_spec S hole).1; omega⟩
    · obtain ⟨x, hx, hxs, hle⟩ := run_taken qs _ _ i q h hq
      exact ⟨x, by simpa using hx, hxs, hle⟩

/-- Later taken slots get strictly larger holes, so no hole is handed out twice. -/
theorem run_strict : ∀ (qs : List ℕ) (hole last i j q q' x y : ℕ), i < j → qs[i]? = some q → q ∈ S →
    qs[j]? = some q' → q' ∈ S → (run S qs hole last)[i]? = some x → (run S qs hole last)[j]? = some y → x < y
  | [], _, _, _, _, _, _, _, _, _, h, _, _, _, _, _ => by simp at h
  | _ :: _, _, _, _, 0, _, _, _, _, hij, _, _, _, _, _, _ => absurd hij (Nat.not_lt_zero _)
  | q0 :: qs, hole, last, 0, j + 1, q, q', x, y, _, hi, hq, hj, hq', hx, hy => by
    simp only [List.getElem?_cons_zero, Option.some.injEq] at hi
    subst hi
    simp only [List.getElem?_cons_succ] at hj
    simp only [run, if_pos hq, List.getElem?_cons_zero, List.getElem?_cons_succ, Option.some.injEq] at hx hy
    subst hx
    obtain ⟨y', hy', -, hle⟩ := run_taken S qs (firstFree S hole + 1) (firstFree S hole) j q' hj hq'
    rw [hy, Option.some.injEq] at hy'
    omega
  | q0 :: qs, hole, last, i + 1, j + 1, q, q', x, y, hij, hi, hq, hj, hq', hx, hy => by
    simp only [List.getElem?_cons_succ] at hi hj
    unfold run at hx hy
    split_ifs at hx hy
    · simp only [List.getElem?_cons_succ] at hx hy
      exact run_strict qs _ _ i j q q' x y (by omega) hi hq hj hq' hx hy
    · simp only [List.getElem?_cons_succ] at hx hy
      exact run_strict qs _ _ i j q q' x y (by omega) hi hq hj hq' hx hy

/-! ## The whole table -/

variable {S}

/-- With n slots in use, as many sit at or past n as there are holes below n. -/
theorem high_eq_free {n T : ℕ} (hS : ∀ q ∈ S, q < T) (hcard : S.card = n) :
    (List.range' n (T - n)).countP (· ∈ S) = free S n := by
  have hhigh : (List.range' n (T - n)).countP (· ∈ S) = (S.filter (n ≤ ·)).card := by
    rw [List.countP_eq_length_filter, ← List.toFinset_card_of_nodup ((List.nodup_range' (s := n) (n := T - n) 1).filter _)]
    congr 1
    ext x
    simp only [List.mem_toFinset, List.mem_filter, List.mem_range'_1, decide_eq_true_eq, Finset.mem_filter]
    constructor
    · rintro ⟨⟨h1, -⟩, h2⟩; exact ⟨h2, h1⟩
    · rintro ⟨h2, h1⟩; exact ⟨⟨h1, by have := hS x h2; omega⟩, h2⟩
  have hlow : (Finset.range n).filter (· ∈ S) = S.filter (¬ n ≤ ·) := by
    ext x; simp only [Finset.mem_filter, Finset.mem_range, not_le]; tauto
  have hsplit := Finset.card_filter_add_card_filter_not (s := S) (n ≤ ·)
  have hrange := Finset.card_filter_add_card_filter_not (s := Finset.range n) (· ∈ S)
  rw [hlow, Finset.card_range] at hrange
  unfold free
  rw [Nat.count_eq_card_filter_range, hhigh]
  omega

/-- Every remap entry is below n, so the hole search never passes n and `EfSeq`'s universe n bounds the table. -/
theorem remap_lt {n T : ℕ} (hS : ∀ q ∈ S, q < T) (hcard : S.card = n) (hn : 0 < n) : ∀ y ∈ remap S n T, y < n :=
  run_lt S n _ 0 0 hn (by rw [high_eq_free hS hcard]; simp [free])

/-- A slot in use at or past n is sent to a free slot below n. -/
theorem lookup_high {n T q : ℕ} (hS : ∀ q ∈ S, q < T) (hcard : S.card = n) (hq : q ∈ S) (hnq : n ≤ q) :
    (remap S n T)[q - n]? = some (lookup S n T q) ∧ lookup S n T q ∉ S ∧ lookup S n T q < n := by
  have hT := hS q hq
  have hidx : (List.range' n (T - n))[q - n]? = some q := by
    rw [List.getElem?_range' (by omega)]; congr 1; omega
  obtain ⟨x, hx, hxs, -⟩ := run_taken S _ 0 0 _ q hidx hq
  have hl : lookup S n T q = x := by simp [lookup, show ¬ q < n by omega, remap, List.getD_eq_getElem?_getD, hx]
  rw [hl]
  exact ⟨hx, hxs, remap_lt hS hcard (hcard ▸ Finset.card_pos.mpr ⟨q, hq⟩) x (List.mem_of_getElem? hx)⟩

theorem lookup_lt {n T q : ℕ} (hS : ∀ q ∈ S, q < T) (hcard : S.card = n) (hq : q ∈ S) : lookup S n T q < n := by
  by_cases h : q < n
  · simp [lookup, h]
  · exact (lookup_high hS hcard hq (by omega)).2.2

theorem lookup_injOn {n T : ℕ} (hS : ∀ q ∈ S, q < T) (hcard : S.card = n) : Set.InjOn (lookup S n T) S := by
  intro q hq q' hq' he
  have hq0 : q ∈ S := hq
  have hq0' : q' ∈ S := hq'
  by_cases h : q < n <;> by_cases h' : q' < n
  · simpa [lookup, h, h'] using he
  · have := (lookup_high hS hcard hq0' (by omega)).2.1
    rw [← he] at this; simp [lookup, h] at this; exact absurd hq0 this
  · have := (lookup_high hS hcard hq0 (by omega)).2.1
    rw [he] at this; simp [lookup, h'] at this; exact absurd hq0' this
  · obtain ⟨hx, -⟩ := lookup_high hS hcard hq0 (T := T) (by omega)
    obtain ⟨hy, -⟩ := lookup_high hS hcard hq0' (T := T) (by omega)
    have hi : (List.range' n (T - n))[q - n]? = some q := by
      rw [List.getElem?_range' (by have := hS q hq0; omega)]; congr 1; omega
    have hj : (List.range' n (T - n))[q' - n]? = some q' := by
      rw [List.getElem?_range' (by have := hS q' hq0'; omega)]; congr 1; omega
    rcases lt_trichotomy (q - n) (q' - n) with hlt | heq | hgt
    · have := run_strict S _ 0 0 _ _ q q' _ _ hlt hi hq0 hj hq0' hx hy; omega
    · omega
    · have := run_strict S _ 0 0 _ _ q' q _ _ hgt hj hq0' hi hq0 hy hx; omega

/-- A successful build is a minimal perfect hash, mapping the n keys one to one onto [0, n).
`pos` is the slot the pilot search placed each key on, injective by `fits_nodup`. -/
theorem lookup_bijective {K : Finset ℕ} {T : ℕ} {pos : ℕ → ℕ} (hinj : Set.InjOn pos K) (hlt : ∀ k ∈ K, pos k < T) :
    Set.InjOn (fun k => lookup (K.image pos) K.card T (pos k)) K ∧
      K.image (fun k => lookup (K.image pos) K.card T (pos k)) = Finset.range K.card := by
  have hS : ∀ q ∈ K.image pos, q < T := by
    intro q hq; obtain ⟨k, hk, rfl⟩ := Finset.mem_image.mp hq; exact hlt k hk
  have hcard : (K.image pos).card = K.card := Finset.card_image_of_injOn hinj
  have hmem : ∀ k ∈ K, pos k ∈ K.image pos := fun k hk => Finset.mem_image_of_mem pos hk
  have hinj' : Set.InjOn (fun k => lookup (K.image pos) K.card T (pos k)) K :=
    fun a ha b hb he => hinj ha hb (lookup_injOn hS hcard (hmem a ha) (hmem b hb) he)
  refine ⟨hinj', Finset.eq_of_subset_of_card_le ?_ ?_⟩
  · intro y hy
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.mp hy
    exact Finset.mem_range.mpr (lookup_lt hS hcard (hmem k hk))
  · rw [Finset.card_image_of_injOn hinj', Finset.card_range]

/-! ## Partitions -/

/-- `part_offset_`, the summed sizes of the partitions before `p`. -/
def offset (sizes : List ℕ) (p : ℕ) : ℕ := (sizes.take p).sum

theorem offset_step {sizes : List ℕ} {p p' : ℕ} (hp : p < sizes.length) (hpp : p < p') :
    offset sizes p + sizes[p] ≤ offset sizes p' := by
  have h1 : offset sizes (p + 1) = offset sizes p + sizes[p] := by
    simp [offset, List.sum_take_succ _ _ hp]
  have h2 : offset sizes p' = offset sizes (p + 1) + ((sizes.drop (p + 1)).take (p' - (p + 1))).sum := by
    rw [offset, offset, show p' = p + 1 + (p' - (p + 1)) by omega, List.take_add, List.sum_append]
    simp
  omega

/-- In `Mphf::operator()`, values from different partitions or different values within one never collide. -/
theorem global_injective {sizes : List ℕ} {p p' l l' : ℕ} (hp : p < sizes.length) (hp' : p' < sizes.length)
    (hl : l < sizes[p]) (hl' : l' < sizes[p']) (he : offset sizes p + l = offset sizes p' + l') : p = p' ∧ l = l' := by
  rcases lt_trichotomy p p' with h | rfl | h
  · have := offset_step hp h; omega
  · exact ⟨rfl, by omega⟩
  · have := offset_step hp' h; omega

/-- Every value `Mphf::operator()` returns is below the total key count. -/
theorem global_lt {sizes : List ℕ} {p l : ℕ} (hp : p < sizes.length) (hl : l < sizes[p]) :
    offset sizes p + l < sizes.sum := by
  have := offset_step hp (show p < p + 1 by omega)
  have h : offset sizes (p + 1) ≤ sizes.sum := by
    rw [offset, ← List.take_append_drop (p + 1) sizes, List.sum_append, List.take_append_drop]; omega
  omega

end OptikProofs.Mphf
