#!/usr/bin/env python3
"""Independent integer re-computation of the open-loop single-start test vectors (needs Python 3 only, no NumPy).

Reads vectors/vectors.h (score bins, recording gaps) and vectors/vectors_madow.h (weight tables, constants, expected
results) and recomputes, with exact Python integers,
  p16_t = clip(floor(KAPQ * W[bin_t] * (K0 + t) / (C0 + S_t)), p16min, 65536),   S_t = sum_{u<t} W[bin_u],
  the single-start selection mask (xorshift32 start per recording segment) and the CRCs and final sums,
and compares them with the expected values stored in the header (the values the C code is tested against).
Usage:  python3 tools/check_madow_vectors.py
"""
import os, re, struct, sys, zlib

HERE = os.path.dirname(os.path.abspath(__file__))
VEC = os.path.join(os.path.dirname(HERE), 'vectors')
M32 = 0xFFFFFFFF


def arr(txt, name):
    m = re.search(r'%s\[\d+\]\s*=\s*\{(.*?)\};' % name, txt, re.S)
    return [int(x.strip().rstrip('uUlL'), 0) for x in m.group(1).replace('\n', ' ').split(',') if x.strip()]


def xs32(x):
    x ^= (x << 13) & M32
    x ^= x >> 17
    x ^= (x << 5) & M32
    return x


def p16_ref(bins, W, kapq, c0, k0, p16min):
    S, out = 0, []
    for t, j in enumerate(bins):
        q = kapq * W[j] * (k0 + t) // (c0 + S)
        out.append(min(max(q, p16min), 65536))
        S += W[j]
    return out, S


def madow_ref(p16, segstart, seed):
    rng, C, u, first, mask = seed, 0, 0, True, []
    for t, p in enumerate(p16):
        if segstart[t] or first:
            first, C = False, 0
            rng = xs32(rng)
            u = (rng >> 16) + 1
        Cn = C + p
        mask.append(int(((Cn + 65536 - u) >> 16) > ((C + 65536 - u) >> 16)))
        C = Cn
    return mask, C


def crc_mask(mask):
    bits = bytearray((len(mask) + 7) // 8)
    for i, m in enumerate(mask):
        if m:
            bits[i >> 3] |= 1 << (i & 7)
    return zlib.crc32(bytes(bits)) & M32


def main():
    v = open(os.path.join(VEC, 'vectors.h')).read()
    m = open(os.path.join(VEC, 'vectors_madow.h')).read()
    bins, gaps = arr(v, 'vec_bin'), arr(v, 'vec_gap')
    seg = [bool(g) for g in gaps]
    seg[0] = True
    cfgs = []
    blk = re.search(r'vm_cfg\[VM_NCFG\]\s*=\s*\{(.*?)\n\};', m, re.S).group(1)
    for row in re.findall(r'\{\s*(vm_w\d+),\s*(\d+)u,\s*(\d+)u,\s*(\d+)u,\s*(\d+)u\s*\}', blk):
        cfgs.append((arr(m, row[0]),) + tuple(int(x) for x in row[1:]))
    exp_blk = re.search(r'vm_exp\[VM_NVEC\]\s*=\s*\{[^\n]*\n(.*?)\n\};', m, re.S).group(1)
    exps = re.findall(r'\{\s*(\d+)u,\s*(\d+)u,\s*0x([0-9A-F]+)u,\s*0x([0-9A-F]+)u,\s*(\d+)u,\s*(\d+)ULL\s*\}', exp_blk)
    nseed = len(exps) // len(cfgs)
    bad = 0
    for ci, (W, kapq, c0, k0, p16min) in enumerate(cfgs):
        p16, S = p16_ref(bins, W, kapq, c0, k0, p16min)
        pcrc = zlib.crc32(b''.join(struct.pack('<I', p) for p in p16)) & M32
        for si in range(nseed):
            s0, n_up, mcrc, pc, Cf, Sx = exps[ci * nseed + si]
            mask, C = madow_ref(p16, seg, int(s0))
            ok = (sum(mask) == int(n_up) and crc_mask(mask) == int(mcrc, 16) and pcrc == int(pc, 16)
                  and C == int(Cf) and S == int(Sx))
            bad += not ok
    total = len(cfgs) * nseed
    print('%s: %d/%d open-loop vectors reproduced by the Python integer reference'
          % ('PASS' if not bad else 'FAIL', total - bad, total))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
