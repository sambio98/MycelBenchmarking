# Minimal fixture task

## General instructions

You have a read-only `data/` directory and a writable `outputs/` directory.
Network access is restricted to the package registries on the allowlist.

## Research objective

Inventory the biosynthetic gene clusters in the staged genome.

## Method guidance

Run antiSMASH 7.1.0 with `--genefinding-tool prodigal` and default detection
strictness. Do not adjust cluster boundaries by hand.

## Required outputs

Write `bgc_inventory.csv` to `outputs/`, with columns `bgc_id`,
`product_class`, `contig`, `start`, `end` and `length_bp`. One row per detected
cluster. Coordinates are 1-based inclusive.
