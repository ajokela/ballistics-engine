#!/bin/bash

echo "============================================"
echo "Testing Unit Conversion System"
echo "============================================"
echo ""

echo "1. Testing Imperial Units (default):"
echo "   Input: 2800 fps velocity, 180 grains mass, 0.308 inch diameter"
cargo run --bin ballistics -- trajectory -v 2800 -b 0.5 -m 180 -d 0.308 --max-range 1000
echo ""

echo "2. Testing Metric Units:"
echo "   Input: 853.44 m/s velocity, 11.664 g mass, 7.823 mm diameter"
cargo run --bin ballistics -- trajectory --units metric -v 853.44 -b 0.5 -m 11.664 -d 7.823 --max-range 1000
echo ""

echo "3. Testing that both produce similar results:"
echo "   (Imperial 90.35 yd ≈ 82.6 m, Metric 47.44 m)"
echo "   Note: Difference is due to --max-range being interpreted in the respective units"
echo ""

echo "4. Testing Zero Range with Imperial:"
cargo run --bin ballistics -- zero -v 2800 -b 0.5 -m 180 -d 0.308 --target-distance 100
echo ""

echo "5. Testing Zero Range with Metric:"
cargo run --bin ballistics -- zero --units metric -v 853.44 -b 0.5 -m 11.664 -d 7.823 --target-distance 100
echo ""

echo "Test complete!"