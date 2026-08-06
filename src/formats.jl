# formats.jl — the 120 exported aliases (ImplementationPlan §5.1).
# GENERATED literal file (plain loop in the source generator, not eval at
# load). Naming per draft §3.2: Binary{K}p{P}{s|u}{e|f}.

"""`Binary3p1sf`: `FloatByte{3,1,SignedFormat,FiniteDomain}` — bitwidth 3, precision 1, signed, finite; exponent bias 2."""
const Binary3p1sf = FloatByte{3,1,SignedFormat,FiniteDomain}
"""`Binary3p1se`: `FloatByte{3,1,SignedFormat,ExtendedDomain}` — bitwidth 3, precision 1, signed, extended (has infinities); exponent bias 2."""
const Binary3p1se = FloatByte{3,1,SignedFormat,ExtendedDomain}
"""`Binary3p2sf`: `FloatByte{3,2,SignedFormat,FiniteDomain}` — bitwidth 3, precision 2, signed, finite; exponent bias 1."""
const Binary3p2sf = FloatByte{3,2,SignedFormat,FiniteDomain}
"""`Binary3p2se`: `FloatByte{3,2,SignedFormat,ExtendedDomain}` — bitwidth 3, precision 2, signed, extended (has infinities); exponent bias 1."""
const Binary3p2se = FloatByte{3,2,SignedFormat,ExtendedDomain}
"""`Binary3p1uf`: `FloatByte{3,1,UnsignedFormat,FiniteDomain}` — bitwidth 3, precision 1, unsigned, finite; exponent bias 4."""
const Binary3p1uf = FloatByte{3,1,UnsignedFormat,FiniteDomain}
"""`Binary3p1ue`: `FloatByte{3,1,UnsignedFormat,ExtendedDomain}` — bitwidth 3, precision 1, unsigned, extended (has infinities); exponent bias 4."""
const Binary3p1ue = FloatByte{3,1,UnsignedFormat,ExtendedDomain}
"""`Binary3p2uf`: `FloatByte{3,2,UnsignedFormat,FiniteDomain}` — bitwidth 3, precision 2, unsigned, finite; exponent bias 2."""
const Binary3p2uf = FloatByte{3,2,UnsignedFormat,FiniteDomain}
"""`Binary3p2ue`: `FloatByte{3,2,UnsignedFormat,ExtendedDomain}` — bitwidth 3, precision 2, unsigned, extended (has infinities); exponent bias 2."""
const Binary3p2ue = FloatByte{3,2,UnsignedFormat,ExtendedDomain}
"""`Binary3p3uf`: `FloatByte{3,3,UnsignedFormat,FiniteDomain}` — bitwidth 3, precision 3, unsigned, finite; exponent bias 1."""
const Binary3p3uf = FloatByte{3,3,UnsignedFormat,FiniteDomain}
"""`Binary3p3ue`: `FloatByte{3,3,UnsignedFormat,ExtendedDomain}` — bitwidth 3, precision 3, unsigned, extended (has infinities); exponent bias 1."""
const Binary3p3ue = FloatByte{3,3,UnsignedFormat,ExtendedDomain}
"""`Binary4p1sf`: `FloatByte{4,1,SignedFormat,FiniteDomain}` — bitwidth 4, precision 1, signed, finite; exponent bias 4."""
const Binary4p1sf = FloatByte{4,1,SignedFormat,FiniteDomain}
"""`Binary4p1se`: `FloatByte{4,1,SignedFormat,ExtendedDomain}` — bitwidth 4, precision 1, signed, extended (has infinities); exponent bias 4."""
const Binary4p1se = FloatByte{4,1,SignedFormat,ExtendedDomain}
"""`Binary4p2sf`: `FloatByte{4,2,SignedFormat,FiniteDomain}` — bitwidth 4, precision 2, signed, finite; exponent bias 2."""
const Binary4p2sf = FloatByte{4,2,SignedFormat,FiniteDomain}
"""`Binary4p2se`: `FloatByte{4,2,SignedFormat,ExtendedDomain}` — bitwidth 4, precision 2, signed, extended (has infinities); exponent bias 2."""
const Binary4p2se = FloatByte{4,2,SignedFormat,ExtendedDomain}
"""`Binary4p3sf`: `FloatByte{4,3,SignedFormat,FiniteDomain}` — bitwidth 4, precision 3, signed, finite; exponent bias 1."""
const Binary4p3sf = FloatByte{4,3,SignedFormat,FiniteDomain}
"""`Binary4p3se`: `FloatByte{4,3,SignedFormat,ExtendedDomain}` — bitwidth 4, precision 3, signed, extended (has infinities); exponent bias 1."""
const Binary4p3se = FloatByte{4,3,SignedFormat,ExtendedDomain}
"""`Binary4p1uf`: `FloatByte{4,1,UnsignedFormat,FiniteDomain}` — bitwidth 4, precision 1, unsigned, finite; exponent bias 8."""
const Binary4p1uf = FloatByte{4,1,UnsignedFormat,FiniteDomain}
"""`Binary4p1ue`: `FloatByte{4,1,UnsignedFormat,ExtendedDomain}` — bitwidth 4, precision 1, unsigned, extended (has infinities); exponent bias 8."""
const Binary4p1ue = FloatByte{4,1,UnsignedFormat,ExtendedDomain}
"""`Binary4p2uf`: `FloatByte{4,2,UnsignedFormat,FiniteDomain}` — bitwidth 4, precision 2, unsigned, finite; exponent bias 4."""
const Binary4p2uf = FloatByte{4,2,UnsignedFormat,FiniteDomain}
"""`Binary4p2ue`: `FloatByte{4,2,UnsignedFormat,ExtendedDomain}` — bitwidth 4, precision 2, unsigned, extended (has infinities); exponent bias 4."""
const Binary4p2ue = FloatByte{4,2,UnsignedFormat,ExtendedDomain}
"""`Binary4p3uf`: `FloatByte{4,3,UnsignedFormat,FiniteDomain}` — bitwidth 4, precision 3, unsigned, finite; exponent bias 2."""
const Binary4p3uf = FloatByte{4,3,UnsignedFormat,FiniteDomain}
"""`Binary4p3ue`: `FloatByte{4,3,UnsignedFormat,ExtendedDomain}` — bitwidth 4, precision 3, unsigned, extended (has infinities); exponent bias 2."""
const Binary4p3ue = FloatByte{4,3,UnsignedFormat,ExtendedDomain}
"""`Binary4p4uf`: `FloatByte{4,4,UnsignedFormat,FiniteDomain}` — bitwidth 4, precision 4, unsigned, finite; exponent bias 1."""
const Binary4p4uf = FloatByte{4,4,UnsignedFormat,FiniteDomain}
"""`Binary4p4ue`: `FloatByte{4,4,UnsignedFormat,ExtendedDomain}` — bitwidth 4, precision 4, unsigned, extended (has infinities); exponent bias 1."""
const Binary4p4ue = FloatByte{4,4,UnsignedFormat,ExtendedDomain}
"""`Binary5p1sf`: `FloatByte{5,1,SignedFormat,FiniteDomain}` — bitwidth 5, precision 1, signed, finite; exponent bias 8."""
const Binary5p1sf = FloatByte{5,1,SignedFormat,FiniteDomain}
"""`Binary5p1se`: `FloatByte{5,1,SignedFormat,ExtendedDomain}` — bitwidth 5, precision 1, signed, extended (has infinities); exponent bias 8."""
const Binary5p1se = FloatByte{5,1,SignedFormat,ExtendedDomain}
"""`Binary5p2sf`: `FloatByte{5,2,SignedFormat,FiniteDomain}` — bitwidth 5, precision 2, signed, finite; exponent bias 4."""
const Binary5p2sf = FloatByte{5,2,SignedFormat,FiniteDomain}
"""`Binary5p2se`: `FloatByte{5,2,SignedFormat,ExtendedDomain}` — bitwidth 5, precision 2, signed, extended (has infinities); exponent bias 4."""
const Binary5p2se = FloatByte{5,2,SignedFormat,ExtendedDomain}
"""`Binary5p3sf`: `FloatByte{5,3,SignedFormat,FiniteDomain}` — bitwidth 5, precision 3, signed, finite; exponent bias 2."""
const Binary5p3sf = FloatByte{5,3,SignedFormat,FiniteDomain}
"""`Binary5p3se`: `FloatByte{5,3,SignedFormat,ExtendedDomain}` — bitwidth 5, precision 3, signed, extended (has infinities); exponent bias 2."""
const Binary5p3se = FloatByte{5,3,SignedFormat,ExtendedDomain}
"""`Binary5p4sf`: `FloatByte{5,4,SignedFormat,FiniteDomain}` — bitwidth 5, precision 4, signed, finite; exponent bias 1."""
const Binary5p4sf = FloatByte{5,4,SignedFormat,FiniteDomain}
"""`Binary5p4se`: `FloatByte{5,4,SignedFormat,ExtendedDomain}` — bitwidth 5, precision 4, signed, extended (has infinities); exponent bias 1."""
const Binary5p4se = FloatByte{5,4,SignedFormat,ExtendedDomain}
"""`Binary5p1uf`: `FloatByte{5,1,UnsignedFormat,FiniteDomain}` — bitwidth 5, precision 1, unsigned, finite; exponent bias 16."""
const Binary5p1uf = FloatByte{5,1,UnsignedFormat,FiniteDomain}
"""`Binary5p1ue`: `FloatByte{5,1,UnsignedFormat,ExtendedDomain}` — bitwidth 5, precision 1, unsigned, extended (has infinities); exponent bias 16."""
const Binary5p1ue = FloatByte{5,1,UnsignedFormat,ExtendedDomain}
"""`Binary5p2uf`: `FloatByte{5,2,UnsignedFormat,FiniteDomain}` — bitwidth 5, precision 2, unsigned, finite; exponent bias 8."""
const Binary5p2uf = FloatByte{5,2,UnsignedFormat,FiniteDomain}
"""`Binary5p2ue`: `FloatByte{5,2,UnsignedFormat,ExtendedDomain}` — bitwidth 5, precision 2, unsigned, extended (has infinities); exponent bias 8."""
const Binary5p2ue = FloatByte{5,2,UnsignedFormat,ExtendedDomain}
"""`Binary5p3uf`: `FloatByte{5,3,UnsignedFormat,FiniteDomain}` — bitwidth 5, precision 3, unsigned, finite; exponent bias 4."""
const Binary5p3uf = FloatByte{5,3,UnsignedFormat,FiniteDomain}
"""`Binary5p3ue`: `FloatByte{5,3,UnsignedFormat,ExtendedDomain}` — bitwidth 5, precision 3, unsigned, extended (has infinities); exponent bias 4."""
const Binary5p3ue = FloatByte{5,3,UnsignedFormat,ExtendedDomain}
"""`Binary5p4uf`: `FloatByte{5,4,UnsignedFormat,FiniteDomain}` — bitwidth 5, precision 4, unsigned, finite; exponent bias 2."""
const Binary5p4uf = FloatByte{5,4,UnsignedFormat,FiniteDomain}
"""`Binary5p4ue`: `FloatByte{5,4,UnsignedFormat,ExtendedDomain}` — bitwidth 5, precision 4, unsigned, extended (has infinities); exponent bias 2."""
const Binary5p4ue = FloatByte{5,4,UnsignedFormat,ExtendedDomain}
"""`Binary5p5uf`: `FloatByte{5,5,UnsignedFormat,FiniteDomain}` — bitwidth 5, precision 5, unsigned, finite; exponent bias 1."""
const Binary5p5uf = FloatByte{5,5,UnsignedFormat,FiniteDomain}
"""`Binary5p5ue`: `FloatByte{5,5,UnsignedFormat,ExtendedDomain}` — bitwidth 5, precision 5, unsigned, extended (has infinities); exponent bias 1."""
const Binary5p5ue = FloatByte{5,5,UnsignedFormat,ExtendedDomain}
"""`Binary6p1sf`: `FloatByte{6,1,SignedFormat,FiniteDomain}` — bitwidth 6, precision 1, signed, finite; exponent bias 16."""
const Binary6p1sf = FloatByte{6,1,SignedFormat,FiniteDomain}
"""`Binary6p1se`: `FloatByte{6,1,SignedFormat,ExtendedDomain}` — bitwidth 6, precision 1, signed, extended (has infinities); exponent bias 16."""
const Binary6p1se = FloatByte{6,1,SignedFormat,ExtendedDomain}
"""`Binary6p2sf`: `FloatByte{6,2,SignedFormat,FiniteDomain}` — bitwidth 6, precision 2, signed, finite; exponent bias 8."""
const Binary6p2sf = FloatByte{6,2,SignedFormat,FiniteDomain}
"""`Binary6p2se`: `FloatByte{6,2,SignedFormat,ExtendedDomain}` — bitwidth 6, precision 2, signed, extended (has infinities); exponent bias 8."""
const Binary6p2se = FloatByte{6,2,SignedFormat,ExtendedDomain}
"""`Binary6p3sf`: `FloatByte{6,3,SignedFormat,FiniteDomain}` — bitwidth 6, precision 3, signed, finite; exponent bias 4."""
const Binary6p3sf = FloatByte{6,3,SignedFormat,FiniteDomain}
"""`Binary6p3se`: `FloatByte{6,3,SignedFormat,ExtendedDomain}` — bitwidth 6, precision 3, signed, extended (has infinities); exponent bias 4."""
const Binary6p3se = FloatByte{6,3,SignedFormat,ExtendedDomain}
"""`Binary6p4sf`: `FloatByte{6,4,SignedFormat,FiniteDomain}` — bitwidth 6, precision 4, signed, finite; exponent bias 2."""
const Binary6p4sf = FloatByte{6,4,SignedFormat,FiniteDomain}
"""`Binary6p4se`: `FloatByte{6,4,SignedFormat,ExtendedDomain}` — bitwidth 6, precision 4, signed, extended (has infinities); exponent bias 2."""
const Binary6p4se = FloatByte{6,4,SignedFormat,ExtendedDomain}
"""`Binary6p5sf`: `FloatByte{6,5,SignedFormat,FiniteDomain}` — bitwidth 6, precision 5, signed, finite; exponent bias 1."""
const Binary6p5sf = FloatByte{6,5,SignedFormat,FiniteDomain}
"""`Binary6p5se`: `FloatByte{6,5,SignedFormat,ExtendedDomain}` — bitwidth 6, precision 5, signed, extended (has infinities); exponent bias 1."""
const Binary6p5se = FloatByte{6,5,SignedFormat,ExtendedDomain}
"""`Binary6p1uf`: `FloatByte{6,1,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 1, unsigned, finite; exponent bias 32."""
const Binary6p1uf = FloatByte{6,1,UnsignedFormat,FiniteDomain}
"""`Binary6p1ue`: `FloatByte{6,1,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 1, unsigned, extended (has infinities); exponent bias 32."""
const Binary6p1ue = FloatByte{6,1,UnsignedFormat,ExtendedDomain}
"""`Binary6p2uf`: `FloatByte{6,2,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 2, unsigned, finite; exponent bias 16."""
const Binary6p2uf = FloatByte{6,2,UnsignedFormat,FiniteDomain}
"""`Binary6p2ue`: `FloatByte{6,2,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 2, unsigned, extended (has infinities); exponent bias 16."""
const Binary6p2ue = FloatByte{6,2,UnsignedFormat,ExtendedDomain}
"""`Binary6p3uf`: `FloatByte{6,3,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 3, unsigned, finite; exponent bias 8."""
const Binary6p3uf = FloatByte{6,3,UnsignedFormat,FiniteDomain}
"""`Binary6p3ue`: `FloatByte{6,3,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 3, unsigned, extended (has infinities); exponent bias 8."""
const Binary6p3ue = FloatByte{6,3,UnsignedFormat,ExtendedDomain}
"""`Binary6p4uf`: `FloatByte{6,4,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 4, unsigned, finite; exponent bias 4."""
const Binary6p4uf = FloatByte{6,4,UnsignedFormat,FiniteDomain}
"""`Binary6p4ue`: `FloatByte{6,4,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 4, unsigned, extended (has infinities); exponent bias 4."""
const Binary6p4ue = FloatByte{6,4,UnsignedFormat,ExtendedDomain}
"""`Binary6p5uf`: `FloatByte{6,5,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 5, unsigned, finite; exponent bias 2."""
const Binary6p5uf = FloatByte{6,5,UnsignedFormat,FiniteDomain}
"""`Binary6p5ue`: `FloatByte{6,5,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 5, unsigned, extended (has infinities); exponent bias 2."""
const Binary6p5ue = FloatByte{6,5,UnsignedFormat,ExtendedDomain}
"""`Binary6p6uf`: `FloatByte{6,6,UnsignedFormat,FiniteDomain}` — bitwidth 6, precision 6, unsigned, finite; exponent bias 1."""
const Binary6p6uf = FloatByte{6,6,UnsignedFormat,FiniteDomain}
"""`Binary6p6ue`: `FloatByte{6,6,UnsignedFormat,ExtendedDomain}` — bitwidth 6, precision 6, unsigned, extended (has infinities); exponent bias 1."""
const Binary6p6ue = FloatByte{6,6,UnsignedFormat,ExtendedDomain}
"""`Binary7p1sf`: `FloatByte{7,1,SignedFormat,FiniteDomain}` — bitwidth 7, precision 1, signed, finite; exponent bias 32."""
const Binary7p1sf = FloatByte{7,1,SignedFormat,FiniteDomain}
"""`Binary7p1se`: `FloatByte{7,1,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 1, signed, extended (has infinities); exponent bias 32."""
const Binary7p1se = FloatByte{7,1,SignedFormat,ExtendedDomain}
"""`Binary7p2sf`: `FloatByte{7,2,SignedFormat,FiniteDomain}` — bitwidth 7, precision 2, signed, finite; exponent bias 16."""
const Binary7p2sf = FloatByte{7,2,SignedFormat,FiniteDomain}
"""`Binary7p2se`: `FloatByte{7,2,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 2, signed, extended (has infinities); exponent bias 16."""
const Binary7p2se = FloatByte{7,2,SignedFormat,ExtendedDomain}
"""`Binary7p3sf`: `FloatByte{7,3,SignedFormat,FiniteDomain}` — bitwidth 7, precision 3, signed, finite; exponent bias 8."""
const Binary7p3sf = FloatByte{7,3,SignedFormat,FiniteDomain}
"""`Binary7p3se`: `FloatByte{7,3,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 3, signed, extended (has infinities); exponent bias 8."""
const Binary7p3se = FloatByte{7,3,SignedFormat,ExtendedDomain}
"""`Binary7p4sf`: `FloatByte{7,4,SignedFormat,FiniteDomain}` — bitwidth 7, precision 4, signed, finite; exponent bias 4."""
const Binary7p4sf = FloatByte{7,4,SignedFormat,FiniteDomain}
"""`Binary7p4se`: `FloatByte{7,4,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 4, signed, extended (has infinities); exponent bias 4."""
const Binary7p4se = FloatByte{7,4,SignedFormat,ExtendedDomain}
"""`Binary7p5sf`: `FloatByte{7,5,SignedFormat,FiniteDomain}` — bitwidth 7, precision 5, signed, finite; exponent bias 2."""
const Binary7p5sf = FloatByte{7,5,SignedFormat,FiniteDomain}
"""`Binary7p5se`: `FloatByte{7,5,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 5, signed, extended (has infinities); exponent bias 2."""
const Binary7p5se = FloatByte{7,5,SignedFormat,ExtendedDomain}
"""`Binary7p6sf`: `FloatByte{7,6,SignedFormat,FiniteDomain}` — bitwidth 7, precision 6, signed, finite; exponent bias 1."""
const Binary7p6sf = FloatByte{7,6,SignedFormat,FiniteDomain}
"""`Binary7p6se`: `FloatByte{7,6,SignedFormat,ExtendedDomain}` — bitwidth 7, precision 6, signed, extended (has infinities); exponent bias 1."""
const Binary7p6se = FloatByte{7,6,SignedFormat,ExtendedDomain}
"""`Binary7p1uf`: `FloatByte{7,1,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 1, unsigned, finite; exponent bias 64."""
const Binary7p1uf = FloatByte{7,1,UnsignedFormat,FiniteDomain}
"""`Binary7p1ue`: `FloatByte{7,1,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 1, unsigned, extended (has infinities); exponent bias 64."""
const Binary7p1ue = FloatByte{7,1,UnsignedFormat,ExtendedDomain}
"""`Binary7p2uf`: `FloatByte{7,2,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 2, unsigned, finite; exponent bias 32."""
const Binary7p2uf = FloatByte{7,2,UnsignedFormat,FiniteDomain}
"""`Binary7p2ue`: `FloatByte{7,2,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 2, unsigned, extended (has infinities); exponent bias 32."""
const Binary7p2ue = FloatByte{7,2,UnsignedFormat,ExtendedDomain}
"""`Binary7p3uf`: `FloatByte{7,3,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 3, unsigned, finite; exponent bias 16."""
const Binary7p3uf = FloatByte{7,3,UnsignedFormat,FiniteDomain}
"""`Binary7p3ue`: `FloatByte{7,3,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 3, unsigned, extended (has infinities); exponent bias 16."""
const Binary7p3ue = FloatByte{7,3,UnsignedFormat,ExtendedDomain}
"""`Binary7p4uf`: `FloatByte{7,4,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 4, unsigned, finite; exponent bias 8."""
const Binary7p4uf = FloatByte{7,4,UnsignedFormat,FiniteDomain}
"""`Binary7p4ue`: `FloatByte{7,4,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 4, unsigned, extended (has infinities); exponent bias 8."""
const Binary7p4ue = FloatByte{7,4,UnsignedFormat,ExtendedDomain}
"""`Binary7p5uf`: `FloatByte{7,5,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 5, unsigned, finite; exponent bias 4."""
const Binary7p5uf = FloatByte{7,5,UnsignedFormat,FiniteDomain}
"""`Binary7p5ue`: `FloatByte{7,5,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 5, unsigned, extended (has infinities); exponent bias 4."""
const Binary7p5ue = FloatByte{7,5,UnsignedFormat,ExtendedDomain}
"""`Binary7p6uf`: `FloatByte{7,6,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 6, unsigned, finite; exponent bias 2."""
const Binary7p6uf = FloatByte{7,6,UnsignedFormat,FiniteDomain}
"""`Binary7p6ue`: `FloatByte{7,6,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 6, unsigned, extended (has infinities); exponent bias 2."""
const Binary7p6ue = FloatByte{7,6,UnsignedFormat,ExtendedDomain}
"""`Binary7p7uf`: `FloatByte{7,7,UnsignedFormat,FiniteDomain}` — bitwidth 7, precision 7, unsigned, finite; exponent bias 1."""
const Binary7p7uf = FloatByte{7,7,UnsignedFormat,FiniteDomain}
"""`Binary7p7ue`: `FloatByte{7,7,UnsignedFormat,ExtendedDomain}` — bitwidth 7, precision 7, unsigned, extended (has infinities); exponent bias 1."""
const Binary7p7ue = FloatByte{7,7,UnsignedFormat,ExtendedDomain}
"""`Binary8p1sf`: `FloatByte{8,1,SignedFormat,FiniteDomain}` — bitwidth 8, precision 1, signed, finite; exponent bias 64."""
const Binary8p1sf = FloatByte{8,1,SignedFormat,FiniteDomain}
"""`Binary8p1se`: `FloatByte{8,1,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 1, signed, extended (has infinities); exponent bias 64."""
const Binary8p1se = FloatByte{8,1,SignedFormat,ExtendedDomain}
"""`Binary8p2sf`: `FloatByte{8,2,SignedFormat,FiniteDomain}` — bitwidth 8, precision 2, signed, finite; exponent bias 32."""
const Binary8p2sf = FloatByte{8,2,SignedFormat,FiniteDomain}
"""`Binary8p2se`: `FloatByte{8,2,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 2, signed, extended (has infinities); exponent bias 32."""
const Binary8p2se = FloatByte{8,2,SignedFormat,ExtendedDomain}
"""`Binary8p3sf`: `FloatByte{8,3,SignedFormat,FiniteDomain}` — bitwidth 8, precision 3, signed, finite; exponent bias 16."""
const Binary8p3sf = FloatByte{8,3,SignedFormat,FiniteDomain}
"""`Binary8p3se`: `FloatByte{8,3,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 3, signed, extended (has infinities); exponent bias 16."""
const Binary8p3se = FloatByte{8,3,SignedFormat,ExtendedDomain}
"""`Binary8p4sf`: `FloatByte{8,4,SignedFormat,FiniteDomain}` — bitwidth 8, precision 4, signed, finite; exponent bias 8."""
const Binary8p4sf = FloatByte{8,4,SignedFormat,FiniteDomain}
"""`Binary8p4se`: `FloatByte{8,4,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 4, signed, extended (has infinities); exponent bias 8."""
const Binary8p4se = FloatByte{8,4,SignedFormat,ExtendedDomain}
"""`Binary8p5sf`: `FloatByte{8,5,SignedFormat,FiniteDomain}` — bitwidth 8, precision 5, signed, finite; exponent bias 4."""
const Binary8p5sf = FloatByte{8,5,SignedFormat,FiniteDomain}
"""`Binary8p5se`: `FloatByte{8,5,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 5, signed, extended (has infinities); exponent bias 4."""
const Binary8p5se = FloatByte{8,5,SignedFormat,ExtendedDomain}
"""`Binary8p6sf`: `FloatByte{8,6,SignedFormat,FiniteDomain}` — bitwidth 8, precision 6, signed, finite; exponent bias 2."""
const Binary8p6sf = FloatByte{8,6,SignedFormat,FiniteDomain}
"""`Binary8p6se`: `FloatByte{8,6,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 6, signed, extended (has infinities); exponent bias 2."""
const Binary8p6se = FloatByte{8,6,SignedFormat,ExtendedDomain}
"""`Binary8p7sf`: `FloatByte{8,7,SignedFormat,FiniteDomain}` — bitwidth 8, precision 7, signed, finite; exponent bias 1."""
const Binary8p7sf = FloatByte{8,7,SignedFormat,FiniteDomain}
"""`Binary8p7se`: `FloatByte{8,7,SignedFormat,ExtendedDomain}` — bitwidth 8, precision 7, signed, extended (has infinities); exponent bias 1."""
const Binary8p7se = FloatByte{8,7,SignedFormat,ExtendedDomain}
"""`Binary8p1uf`: `FloatByte{8,1,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 1, unsigned, finite; exponent bias 128."""
const Binary8p1uf = FloatByte{8,1,UnsignedFormat,FiniteDomain}
"""`Binary8p1ue`: `FloatByte{8,1,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 1, unsigned, extended (has infinities); exponent bias 128."""
const Binary8p1ue = FloatByte{8,1,UnsignedFormat,ExtendedDomain}
"""`Binary8p2uf`: `FloatByte{8,2,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 2, unsigned, finite; exponent bias 64."""
const Binary8p2uf = FloatByte{8,2,UnsignedFormat,FiniteDomain}
"""`Binary8p2ue`: `FloatByte{8,2,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 2, unsigned, extended (has infinities); exponent bias 64."""
const Binary8p2ue = FloatByte{8,2,UnsignedFormat,ExtendedDomain}
"""`Binary8p3uf`: `FloatByte{8,3,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 3, unsigned, finite; exponent bias 32."""
const Binary8p3uf = FloatByte{8,3,UnsignedFormat,FiniteDomain}
"""`Binary8p3ue`: `FloatByte{8,3,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 3, unsigned, extended (has infinities); exponent bias 32."""
const Binary8p3ue = FloatByte{8,3,UnsignedFormat,ExtendedDomain}
"""`Binary8p4uf`: `FloatByte{8,4,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 4, unsigned, finite; exponent bias 16."""
const Binary8p4uf = FloatByte{8,4,UnsignedFormat,FiniteDomain}
"""`Binary8p4ue`: `FloatByte{8,4,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 4, unsigned, extended (has infinities); exponent bias 16."""
const Binary8p4ue = FloatByte{8,4,UnsignedFormat,ExtendedDomain}
"""`Binary8p5uf`: `FloatByte{8,5,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 5, unsigned, finite; exponent bias 8."""
const Binary8p5uf = FloatByte{8,5,UnsignedFormat,FiniteDomain}
"""`Binary8p5ue`: `FloatByte{8,5,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 5, unsigned, extended (has infinities); exponent bias 8."""
const Binary8p5ue = FloatByte{8,5,UnsignedFormat,ExtendedDomain}
"""`Binary8p6uf`: `FloatByte{8,6,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 6, unsigned, finite; exponent bias 4."""
const Binary8p6uf = FloatByte{8,6,UnsignedFormat,FiniteDomain}
"""`Binary8p6ue`: `FloatByte{8,6,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 6, unsigned, extended (has infinities); exponent bias 4."""
const Binary8p6ue = FloatByte{8,6,UnsignedFormat,ExtendedDomain}
"""`Binary8p7uf`: `FloatByte{8,7,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 7, unsigned, finite; exponent bias 2."""
const Binary8p7uf = FloatByte{8,7,UnsignedFormat,FiniteDomain}
"""`Binary8p7ue`: `FloatByte{8,7,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 7, unsigned, extended (has infinities); exponent bias 2."""
const Binary8p7ue = FloatByte{8,7,UnsignedFormat,ExtendedDomain}
"""`Binary8p8uf`: `FloatByte{8,8,UnsignedFormat,FiniteDomain}` — bitwidth 8, precision 8, unsigned, finite; exponent bias 1."""
const Binary8p8uf = FloatByte{8,8,UnsignedFormat,FiniteDomain}
"""`Binary8p8ue`: `FloatByte{8,8,UnsignedFormat,ExtendedDomain}` — bitwidth 8, precision 8, unsigned, extended (has infinities); exponent bias 1."""
const Binary8p8ue = FloatByte{8,8,UnsignedFormat,ExtendedDomain}

"""
    FORMAT_ALIASES

All 120 format aliases in canonical lattice order (`format_id` order):
K ascending; Signed before Unsigned; P ascending; Finite before Extended.
"""
const FORMAT_ALIASES = (
    Binary3p1sf,
    Binary3p1se,
    Binary3p2sf,
    Binary3p2se,
    Binary3p1uf,
    Binary3p1ue,
    Binary3p2uf,
    Binary3p2ue,
    Binary3p3uf,
    Binary3p3ue,
    Binary4p1sf,
    Binary4p1se,
    Binary4p2sf,
    Binary4p2se,
    Binary4p3sf,
    Binary4p3se,
    Binary4p1uf,
    Binary4p1ue,
    Binary4p2uf,
    Binary4p2ue,
    Binary4p3uf,
    Binary4p3ue,
    Binary4p4uf,
    Binary4p4ue,
    Binary5p1sf,
    Binary5p1se,
    Binary5p2sf,
    Binary5p2se,
    Binary5p3sf,
    Binary5p3se,
    Binary5p4sf,
    Binary5p4se,
    Binary5p1uf,
    Binary5p1ue,
    Binary5p2uf,
    Binary5p2ue,
    Binary5p3uf,
    Binary5p3ue,
    Binary5p4uf,
    Binary5p4ue,
    Binary5p5uf,
    Binary5p5ue,
    Binary6p1sf,
    Binary6p1se,
    Binary6p2sf,
    Binary6p2se,
    Binary6p3sf,
    Binary6p3se,
    Binary6p4sf,
    Binary6p4se,
    Binary6p5sf,
    Binary6p5se,
    Binary6p1uf,
    Binary6p1ue,
    Binary6p2uf,
    Binary6p2ue,
    Binary6p3uf,
    Binary6p3ue,
    Binary6p4uf,
    Binary6p4ue,
    Binary6p5uf,
    Binary6p5ue,
    Binary6p6uf,
    Binary6p6ue,
    Binary7p1sf,
    Binary7p1se,
    Binary7p2sf,
    Binary7p2se,
    Binary7p3sf,
    Binary7p3se,
    Binary7p4sf,
    Binary7p4se,
    Binary7p5sf,
    Binary7p5se,
    Binary7p6sf,
    Binary7p6se,
    Binary7p1uf,
    Binary7p1ue,
    Binary7p2uf,
    Binary7p2ue,
    Binary7p3uf,
    Binary7p3ue,
    Binary7p4uf,
    Binary7p4ue,
    Binary7p5uf,
    Binary7p5ue,
    Binary7p6uf,
    Binary7p6ue,
    Binary7p7uf,
    Binary7p7ue,
    Binary8p1sf,
    Binary8p1se,
    Binary8p2sf,
    Binary8p2se,
    Binary8p3sf,
    Binary8p3se,
    Binary8p4sf,
    Binary8p4se,
    Binary8p5sf,
    Binary8p5se,
    Binary8p6sf,
    Binary8p6se,
    Binary8p7sf,
    Binary8p7se,
    Binary8p1uf,
    Binary8p1ue,
    Binary8p2uf,
    Binary8p2ue,
    Binary8p3uf,
    Binary8p3ue,
    Binary8p4uf,
    Binary8p4ue,
    Binary8p5uf,
    Binary8p5ue,
    Binary8p6uf,
    Binary8p6ue,
    Binary8p7uf,
    Binary8p7ue,
    Binary8p8uf,
    Binary8p8ue,
)

export Binary3p1sf
export Binary3p1se
export Binary3p2sf
export Binary3p2se
export Binary3p1uf
export Binary3p1ue
export Binary3p2uf
export Binary3p2ue
export Binary3p3uf
export Binary3p3ue
export Binary4p1sf
export Binary4p1se
export Binary4p2sf
export Binary4p2se
export Binary4p3sf
export Binary4p3se
export Binary4p1uf
export Binary4p1ue
export Binary4p2uf
export Binary4p2ue
export Binary4p3uf
export Binary4p3ue
export Binary4p4uf
export Binary4p4ue
export Binary5p1sf
export Binary5p1se
export Binary5p2sf
export Binary5p2se
export Binary5p3sf
export Binary5p3se
export Binary5p4sf
export Binary5p4se
export Binary5p1uf
export Binary5p1ue
export Binary5p2uf
export Binary5p2ue
export Binary5p3uf
export Binary5p3ue
export Binary5p4uf
export Binary5p4ue
export Binary5p5uf
export Binary5p5ue
export Binary6p1sf
export Binary6p1se
export Binary6p2sf
export Binary6p2se
export Binary6p3sf
export Binary6p3se
export Binary6p4sf
export Binary6p4se
export Binary6p5sf
export Binary6p5se
export Binary6p1uf
export Binary6p1ue
export Binary6p2uf
export Binary6p2ue
export Binary6p3uf
export Binary6p3ue
export Binary6p4uf
export Binary6p4ue
export Binary6p5uf
export Binary6p5ue
export Binary6p6uf
export Binary6p6ue
export Binary7p1sf
export Binary7p1se
export Binary7p2sf
export Binary7p2se
export Binary7p3sf
export Binary7p3se
export Binary7p4sf
export Binary7p4se
export Binary7p5sf
export Binary7p5se
export Binary7p6sf
export Binary7p6se
export Binary7p1uf
export Binary7p1ue
export Binary7p2uf
export Binary7p2ue
export Binary7p3uf
export Binary7p3ue
export Binary7p4uf
export Binary7p4ue
export Binary7p5uf
export Binary7p5ue
export Binary7p6uf
export Binary7p6ue
export Binary7p7uf
export Binary7p7ue
export Binary8p1sf
export Binary8p1se
export Binary8p2sf
export Binary8p2se
export Binary8p3sf
export Binary8p3se
export Binary8p4sf
export Binary8p4se
export Binary8p5sf
export Binary8p5se
export Binary8p6sf
export Binary8p6se
export Binary8p7sf
export Binary8p7se
export Binary8p1uf
export Binary8p1ue
export Binary8p2uf
export Binary8p2ue
export Binary8p3uf
export Binary8p3ue
export Binary8p4uf
export Binary8p4ue
export Binary8p5uf
export Binary8p5ue
export Binary8p6uf
export Binary8p6ue
export Binary8p7uf
export Binary8p7ue
export Binary8p8uf
export Binary8p8ue
