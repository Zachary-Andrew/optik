# Appendix B: formal statements

Each row ties a Lean theorem to the paper item it backs and the code it models.
Status is "proved" only when `lake build` accepts the theorem and `python3 audit.py` finds no `sorry` and only Lean's standard axioms (`propext`, `Classical.choice`, `Quot.sound`) behind it.
Theorems proved with `bv_decide`, and the theorems that use them, also depend on the `_native.bv_decide` axioms those calls register.
Each records that Lean's compiled checker verified one SAT certificate, because the kernel cannot evaluate 64-bit words directly.
`lean/audit.py` lists the theorems allowed to register them (`NATIVE_SOURCES`) and fails on a native axiom from anywhere else.

Lean checks these statements about the algorithm model.
The last column says what ties the model to the C++ and Rust.
"Vectors" means the model's test vectors are run through the binary by `bench/validate_lean.py` (step 12 of `bench/run_validations.sh`) or by a Rust unit test.
"Direct" means the Lean definition transcribes the code's own operations, so the theorem is about that code path.

| Lean theorem | Paper item | Models | Assumptions | Status | Tie to the code |
|---|---|---|---|---|---|
| `Kmer.comp_val_eq_xor`, `Kmer.rc_rc` | Lemma 4.1 | `base_complement`, `reverse_complement`, `include/kmer.hpp` | none | proved | direct |
| `Kmer.pack_injective`, `Kmer.canon_rc`, `Kmer.canon_eq_iff` | Lemma 4.1 | `encode_kmer` in `include/kmer.hpp`, and the `canon` members (for example `include/dbg_part.hpp`) | k-mers of one length | proved | direct |
| `Kmer.rcWord_eq_rcSpec` | Lemma 4.1 | the bit trick in `reverse_complement`, `include/kmer.hpp` | k ≤ 32, native SAT-certificate check | proved | direct, for every 64-bit input |
| `Hash.fmix64_injective`, `Hash.mhash_ordInj` | discharges `Spss.OrdInj` | `mix64` in `include/wyhash.hpp`, `SpssCore::mhash` | m ≤ 32, native SAT-certificate check | proved | direct |
| `Hash.wang_injOn` | "distinct k-mers can never tie" | `minimizer_order_hash`, `include/minimizer.hpp` | mask of n ≤ 64 bits, native SAT-certificate check | proved | direct |
| `Rolling.run_spec`, `Rolling.run_none` | Lemma 4.1 (extraction) | `extract_kmers`, `include/kmer_extractor.hpp`; the same updates in `src/kmer_count.cpp` | 1 ≤ k ≤ 32, native SAT-certificate check | proved | direct |
| `Sort.parSort_eq_ref` | Lemma 4.3 | `par_sort`, `include/par_sort.hpp` | linear order, any slicing | proved | direct |
| `Sort.radixSort_eq_ref` | Lemma 4.3 | `radix_sort_inplace`, `include/par_sort.hpp` | bucket key monotone in the order | proved | direct |
| `Minimizer.stream_isLeftMin`, `Minimizer.stream_eq_minPos` | Lemma 4.2 (window rule) | `extract_minimizers`, `include/minimizer.hpp` | one N-free run | proved | direct |
| `Minimizer.minPos_monotone` | used by prop:seed | `extract_minimizers` | one N-free run | proved | direct |
| `Minimizer.tied_mirror`, `Minimizer.order_mirror` | Lemma 4.2 | SpssCore's m-mer minimizer | order depends only on the canonical m-mer | proved | through the Spss vectors |
| `SeedIndex.parallel_eq_serial` | prop:seed | step 1 of `SeedIndex::build`, `include/seed_index.hpp` | one N-free run, cut points inside the reference | proved | direct |
| `Spss.present_iff` | prop:spss | lookup, `include/spss_core.hpp` | MPHF injective on its build set (`Mphf.lookup_bijective`, for a build that returned true), order injective on packed m-mers (`Spss.OrdInj`, proved for the code's order by `Hash.mhash_ordInj`), candidates include every verifying record | proved | vectors (`build-spss`, `query --ids`) |
| `Spss.extension_sound` | prop:spss (streaming) | extension, `include/spss_scan.hpp` | as `present_iff` | proved | vectors |
| `Spss.ids_dense` | dense ids | `kmer_id`, `include/spss_core.hpp` | k ≥ 1, every stored string at least k long | proved | direct |
| `Spss.occursB_iff` | vector oracle | the Spss vectors' expected answers | none | proved | — |
| `Dbg.compacts_rc` | Definition 7.1 | the walker's rule, `include/dbg_unitigs.hpp` | symmetric edge masks (`Count.code_symmetric`, `Count.trim_symmetric`) | proved | vectors |
| `Dbg.unitig_starts`, `Dbg.unitig_rc` | Lemma 7.2 | `run_sorted`, `include/dbg_unitigs.hpp` | symmetric edge masks (as `compacts_rc`) | proved | vectors (`--graph table`, and at -c 2 `dbg2.json`) |
| `Dbg.glue_iff_compacts`, `Dbg.partStep_eq_compacts` | prop:glue | `PartBuilder`, `include/dbg_part.hpp` | partition is a function of the canonical k-mer (`Count.minKey_canon`) | proved | vectors (`--graph part` vs `--graph table`) |
| `Dbg.rightKey_eq_leftKey`, `Dbg.key_unique` | prop:glue (keys) | open-end keys, `include/dbg_part.hpp` | k ≥ 1 for `key_unique`, overlapping edges | proved | vectors |
| `Dbg.claimers_are_endpoints` | cutoff claim protocol | cross-partition claims, `include/dbg_part.hpp` | k ≥ 1, overlapping edges | proved | — |
| `Mapq.posterior_identity` | prop:post | `posterior_mapq_w`, `tools/src/map/mapq.rs` | equal priors, likelihood c · r^d | proved | Rust test against vectors |
| `Mapq.mapq_antitone_S`, `Mapq.mapq_antitone_H` | prop:post (monotonicity) | `posterior_mapq_w` | positive unseen-locus odds | proved | Rust test against vectors |
| `Mapq.tie_without_override`, `Mapq.tie_with_override`, `Mapq.phred_half` | audit item 4 | the MAPQ 0 rule in `posterior_mapq_w` | none | proved | Rust test against vectors |
| `Mapq.calibration_bound` | prop:cal | the calibration test's bound | Le Cam (1960) as a named hypothesis | proved | — |
| `Align.shift_deletion_le`, `Align.shift_insertion_le` | Prop. 6, Prop. 6′ | `left_align`, `tools/src/map/cigar.rs` | per-column mismatch cost, unit gap cost (NM) | proved | direct |
| `Dyn.insert_set`, `Dyn.erase_set`, `Dyn.insert_result`, `Dyn.erase_result`, `Dyn.compact_set` | editable index exactness | `DynSpss`, `include/dyn_spss.hpp`; `okdyn_apply_fasta` in `src/ffi_spss.cpp` | edit buffers modelled as finite sets, keys canonical | proved | direct |
| `Varint.decode_encode`, `Varint.decodeNodes_encodeNodes` | spill exactness | `put_varint`, `save_nodes` and the reader in `src/kmer_count.cpp` | native SAT-certificate check | proved | direct, and vectors (`dbg2.json` with `OPTIK_SPILL_TEST`) |
| `Paired.frTlen_comm`, `Paired.proper_comm` | mate-order independence | `fr_tlen`, `proper`, `tools/src/map/paired.rs` | none | proved | direct |
| `Paired.tlen_opposite`, `Paired.tlen_left_positive`, `Paired.tlen_ne_zero_iff` | SAM TLEN rules | `set_pair_fields`, `tools/src/map/paired.rs` | each record's first base at or before its last | proved | direct |
| `Mphf.fits_nodup`, `Mphf.lookup_bijective` | MPHF is a bijection onto [0, n) | pilot search, remap and `operator()`, `PilotMphf` in `include/phobic_mphf.hpp` | a build that returned true | proved | direct |
| `Mphf.remap_sorted`, `Mphf.remap_lt` | the remap table suits `EfSeq`, and the hole search never passes n | remap loop in `PilotMphf::try_build` | as `lookup_bijective` | proved | direct |
| `Mphf.global_injective`, `Mphf.global_lt` | partitioned values are distinct and below n | `Mphf::operator()`, `part_offset_` | none | proved | direct |
| `Cigar.expand_pushOp`, `Cigar.canon_pushOp`, `Cigar.pushAll_spec` | merged CIGARs have no empty runs or equal neighbours | `push_op`, `push_cigar`, `tools/src/map/cigar.rs` | none | proved | direct |
| `Cigar.refSpan_eq`, `Cigar.colCount_eq`, `Cigar.qlen_eq` | END, TLEN and column counts | `ref_span` (`paired.rs`), `cigar_cols` | none | proved | vectors (`cigar.json`, Rust `cigar_helpers_match_lean_vectors`) |
| `Cigar.edImplied_eq_cost` | NM is the CIGAR's edit cost | `cigar_ed_implied` | none | proved | vectors (`cigar.json`) |
| `Cigar.clip_query_length`, `Cigar.clip_nm`, `Cigar.clip_nm_le` | soft-clipped records have consistent SEQ and NM | `clip_ends` | any cut points with `lcut + rcut` below the column count, as the code requires | proved | direct |
| `Cigar.hardClip_qlen`, `Cigar.hardClip_refSpan` | hard clips leave SEQ length and POS unchanged | `hard_clip` | none (an absent clip is a run of length 0) | proved | direct |
| `Count.slide_min`, `Count.slide_eq_minKey` | super-k-mer minimizers | the rolling window minimum in pass 1 of `count_into_graph`, `src/kmer_count.cpp` | order values fit 64 bits | proved | direct |
| `Count.minOrd_rc`, `Count.minKey_canon` | partition is strand-blind | `sk_min` and the partition hash | m ≤ k, k-mers of length k | proved | direct |
| `Count.records_flatten`, `Count.records_key`, `Count.count_partition` | partitioned counting is exact | the record emit loop and pass 2's per-partition tables | the key depends only on the canonical k-mer (`minKey_canon`) | proved | vectors (`dbg2.json`, -c 2) |
| `Count.bits_iff_obs`, `Count.codeEdge_iff` | edge masks are the observed (k+1)-mers | `in_b`/`out_b` in pass 2's `each_kmer` | each occurrence is read with its read neighbours (the flank codes), node bits in canonical orientation | proved | vectors (`dbg2.json`, GFA edges) |
| `Count.obs_symmetric`, `Count.code_symmetric`, `Count.trim_symmetric` | discharges the Dbg `hsym` | the masks and the solid-neighbour trim | k ≥ 1 | proved | vectors (`dbg2.json`, GFA edges) |
| `Bpm.row`, `Bpm.bounds`, `Bpm.cell` | the delta form of the edit DP | Myers' recurrence, `include/bpm_edit.hpp` | costs 0 or 1 | proved | direct |
| `Bpm.carry`, `Bpm.wordStep_bits` | the word step computes 64 rows | the loop body of `bpm_semiglobal_ed` | native SAT-certificate check | proved | direct |
| `Bpm.column`, `Bpm.run_eq_dp`, `Bpm.bpm_eq_dp` | placement verification is exact | `bpm_semiglobal_ed`, several words chained through `hin` | 0 < m, the code's word count ⌈m / 64⌉, match words as the code builds them | proved | vectors (`bpm.json`, Rust `placements_match_lean_vectors`, scalar and AVX2 lanes) |
| `Ef.key_eq` | node i is the i-th smallest k-mer | `DbgGraphEF::key`, `EfSeq::operator[]` | keys in non-decreasing order, select1 as specified | proved | direct |
| `Ef.sel0_bucket`, `Ef.rank_some`, `Ef.rank_complete` | Elias-Fano membership is exact | `DbgGraphEF::rank`, `include/dbg_graph_ef.hpp` | keys strictly increasing, select0 as specified (the sampled word walk is not modelled) | proved | vectors (`dbg.json`, `dbg2.json` with `OPTIK_DBG_COMPACT=1`) |
| `Fm.occs_sep`, `Fm.count_text`, `Fm.occurs_text` | prop:kbwt (the windows half) | the X-separated text of `FmIndex`, `include/fm_index.hpp` | the query is ACGT only and non-empty | proved | direct |
| `FmIndex.*` | prop:kbwt | `include/fm_index.hpp` | — | stretch goal, not started | — |

Helper lemmas in the same files (for example `Sort.mergeAll_spec` or `Spss.key_rc`) are proved under the same audit and serve the rows above.
