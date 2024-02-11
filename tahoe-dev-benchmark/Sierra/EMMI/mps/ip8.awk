BEGIN {skip = 0}
skip == 7 { print $4,$5,$6,$7,$8,$9,$10,$11,$12; skip = 0; next }
skip < 7 { skip++; next }
