#!/bin/bash
echo "Testing imperial units (default):"
cargo run --bin ballistics -- trajectory \
  -v 2800 \
  -b 0.5 \
  -m 180 \
  -d 0.308 \
  --max-range 1000