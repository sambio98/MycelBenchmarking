"""Benchmark the Biomni biomedical agent over a set of Q&A pairs.

A single CLI command (``biomni-benchmark``) runs every question in a JSONL file
through Biomni and collects the raw outputs for later scoring. The underlying
LLM is fully configurable so Biomni can be tested on different models.
"""

__version__ = "0.1.0"
