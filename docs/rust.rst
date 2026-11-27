Rust front end
==============

The ``optik`` binary parses arguments, drives the C++ core through the C ABI, writes SAM records and the metric lines, and hosts the interactive TUI.
Its API reference is rustdoc's output, private items included, since nothing in the crate is a public library.

`Open the Rust API reference <rust/optik/index.html>`_

``docs/build.sh`` builds it with ``cargo doc --no-deps --document-private-items`` and copies ``target/doc`` into the site.
