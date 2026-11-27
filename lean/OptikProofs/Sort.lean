import Mathlib.Data.List.Sort
import Mathlib.Data.List.FinRange
import Mathlib.Algebra.BigOperators.Fin

/-!
# Deterministic parallel sorting

Models `include/par_sort.hpp`.
Both sorts return exactly what `std::sort` returns under a strict total order, which makes the index builders' output identical at every thread count.
The argument is the same for both.
Each produces a sorted permutation of its input, and under a linear order a list has only one sorted permutation.
`ref`, a plain merge sort, stands for `std::sort`.
-/

namespace OptikProofs.Sort

variable {α : Type*} [LinearOrder α]

/-- The reference result that every parallel sort must reproduce. -/
def ref (l : List α) : List α := l.mergeSort (fun a b => a ≤ b)

theorem ref_perm (l : List α) : (ref l).Perm l := List.mergeSort_perm l _

theorem ref_sorted (l : List α) : (ref l).Pairwise (· ≤ ·) := List.pairwise_mergeSort' _ l

/-- Under a linear order a list has one sorted permutation, so any sort that outputs a sorted permutation returns `ref`. -/
theorem eq_ref {l out : List α} (hp : out.Perm l) (hs : out.Pairwise (· ≤ ·)) : out = ref l :=
  List.Perm.eq_of_pairwise' hs (ref_sorted l) (hp.trans (ref_perm l).symm)

/-! ## `par_sort` merges sorted slices pairwise -/

/-- One round of `par_sort`'s merge loop merges neighbouring runs in pairs, and an odd last run passes to the next round unchanged. -/
def mergeRound : List (List α) → List (List α)
  | a :: b :: rest => List.merge a b (fun x y => x ≤ y) :: mergeRound rest
  | runs => runs

theorem mergeRound_length_lt : ∀ (a b : List α) (rest : List (List α)),
    (mergeRound (a :: b :: rest)).length < (a :: b :: rest).length
  | a, b, [] => by simp [mergeRound]
  | a, b, [c] => by simp [mergeRound]
  | a, b, c :: d :: rest => by
    have := mergeRound_length_lt c d rest
    simp only [mergeRound, List.length_cons] at this ⊢
    omega

/-- Rounds repeat until at most one run is left, as the width-doubling loop in `par_sort` does. -/
def mergeAll : List (List α) → List α
  | [] => []
  | [run] => run
  | a :: b :: rest =>
    have := mergeRound_length_lt a b rest
    mergeAll (mergeRound (a :: b :: rest))
termination_by runs => runs.length

/-- `par_sort` on input split into contiguous slices `slices`, each sorted by one thread before the merge rounds. -/
def parSort (slices : List (List α)) : List α := mergeAll (slices.map ref)

/-- Each run handed to the merge rounds is sorted. -/
def AllSorted (runs : List (List α)) : Prop := ∀ r ∈ runs, r.Pairwise (· ≤ ·)

theorem mergeRound_sorted : ∀ runs : List (List α), AllSorted runs → AllSorted (mergeRound runs)
  | [], h => h
  | [_], h => h
  | a :: b :: rest, h => by
    intro r hr
    simp only [mergeRound, List.mem_cons] at hr
    rcases hr with rfl | hr
    · exact (h a (by simp)).merge (h b (by simp))
    · exact mergeRound_sorted rest (fun r' hr' => h r' (by simp [hr'])) r hr

/-- A merge round moves elements between runs but neither adds nor drops any. -/
theorem mergeRound_perm : ∀ runs : List (List α), (mergeRound runs).flatten.Perm runs.flatten
  | [] => by simp [mergeRound]
  | [_] => by simp [mergeRound]
  | a :: b :: rest => by
    simp only [mergeRound, List.flatten_cons, ← List.append_assoc]
    exact List.Perm.append (List.merge_perm_append _) (mergeRound_perm rest)

theorem mergeAll_spec (runs : List (List α)) (h : AllSorted runs) :
    (mergeAll runs).Pairwise (· ≤ ·) ∧ (mergeAll runs).Perm runs.flatten := by
  induction runs using mergeAll.induct with
  | case1 => simp [mergeAll]
  | case2 run => exact ⟨by simpa [mergeAll] using h run (by simp), by simp [mergeAll]⟩
  | case3 a b rest _ ih =>
    obtain ⟨hs, hp⟩ := ih (mergeRound_sorted _ h)
    rw [mergeAll]
    exact ⟨hs, hp.trans (mergeRound_perm _)⟩

/-- `par_sort` returns what `std::sort` returns, whatever slices the threads were given. -/
theorem parSort_eq_ref (slices : List (List α)) : parSort slices = ref slices.flatten := by
  have hperm : (slices.map ref).flatten.Perm slices.flatten := by
    induction slices with
    | nil => simp
    | cons s rest ih => simpa using (ref_perm s).append ih
  have hsorted : AllSorted (slices.map ref) := by
    intro r hr
    obtain ⟨s, -, rfl⟩ := List.mem_map.mp hr
    exact ref_sorted s
  obtain ⟨hs, hp⟩ := mergeAll_spec _ hsorted
  exact eq_ref (hp.trans hperm) hs

/-! ## `radix_sort_inplace` distributes once, then sorts each bucket -/

section Radix

variable {B : ℕ} (key : α → Fin B)

/-- After the American-flag pass, bucket j holds, in input order, the elements whose top key bits equal j. -/
def buckets (l : List α) : List (List α) := (List.finRange B).map fun j => l.filter fun a => key a = j

/-- `radix_sort_inplace` sorts every bucket on its own and lays the buckets out in bucket order. -/
def radixSort (l : List α) : List α := ((buckets key l).map ref).flatten

/-- An element shows up only in the bucket of its own key. -/
theorem count_filter_key (l : List α) (a : α) (j : Fin B) :
    (l.filter fun b => decide (key b = j)).count a = if key a = j then l.count a else 0 := by
  split_ifs with h
  · exact List.count_filter (by simpa using h)
  · exact List.count_eq_zero.mpr (by simp [List.mem_filter, h])

/-- The buckets together hold every element exactly as often as the input does. -/
theorem flatten_buckets_perm (l : List α) : ((buckets key l).map ref).flatten.Perm l := by
  rw [List.perm_iff_count]
  intro a
  simp only [buckets, List.map_map, List.count_flatten, Function.comp_def, (ref_perm _).count_eq,
    count_filter_key]
  rw [← Fin.sum_univ_def]
  simp

/-- `radix_sort_inplace` returns what `std::sort` returns, provided the comparator orders by the bucket key first, as par_sort.hpp requires. -/
theorem radixSort_eq_ref (hkey : Monotone key) (l : List α) : radixSort key l = ref l := by
  refine eq_ref (flatten_buckets_perm key l) ?_
  rw [radixSort, List.pairwise_flatten]
  refine ⟨?_, ?_⟩
  · intro r hr
    obtain ⟨s, -, rfl⟩ := List.mem_map.mp hr
    exact ref_sorted s
  · simp only [buckets, List.map_map, List.pairwise_map, Function.comp_def]
    refine (List.pairwise_lt_finRange B).imp ?_
    intro i j hij x hx y hy
    have hxk : key x = i := (by simpa using (ref_perm _).mem_iff.mp hx : x ∈ l ∧ key x = i).2
    have hyk : key y = j := (by simpa using (ref_perm _).mem_iff.mp hy : y ∈ l ∧ key y = j).2
    by_contra hxy
    have := hkey (le_of_lt (not_le.mp hxy))
    rw [hxk, hyk] at this
    exact absurd hij (not_lt.mpr this)

end Radix

end OptikProofs.Sort
