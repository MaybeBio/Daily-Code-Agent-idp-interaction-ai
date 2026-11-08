"""
xtc_stdlib.py — read GROMACS/LAMMPS .xtc trajectories with the Python standard library only
(mdtraj cannot be built on this Mac: the Command Line Tools linker rejects the macOS SDK).

A port of the decompression in xdrfile.c (xdr3dfcoord): big-endian XDR header, then integer coordinates packed
with the "magicints" scheme (runs of small differences between neighbouring atoms). Coordinates are returned in nm.

  frames = index(data)                  # byte offset of every frame (cheap: reads headers and block sizes only)
  natoms, step, time, box, X, Y, Z = decode(data, frames[k])
"""
import struct

MAGICINTS = [0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
             80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
             1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
             16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
             131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
             832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
             4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216]
FIRSTIDX = 9
MAGIC = 1995


def _header(data, off):
    magic, natoms, step = struct.unpack_from('>iii', data, off)
    if magic != MAGIC:
        raise ValueError(f'bad xtc magic {magic} at byte {off}')
    time = struct.unpack_from('>f', data, off + 12)[0]
    box = struct.unpack_from('>9f', data, off + 16)
    lsize = struct.unpack_from('>i', data, off + 52)[0]
    return natoms, step, time, box, lsize, off + 56


def index(data):
    """Byte offsets of all complete frames."""
    offs = []
    off, n = 0, len(data)
    while off + 56 <= n:
        natoms, step, time, box, lsize, p = _header(data, off)
        if lsize <= 9:
            end = p + 12 * lsize
        else:
            nbytes = struct.unpack_from('>i', data, p + 32)[0]
            end = p + 36 + ((nbytes + 3) // 4) * 4
        if end > n:
            break
        offs.append(off)
        off = end
    return offs


class _Bits:
    """MSB-first bit reader over a byte string (receivebits in xdrfile.c)."""
    __slots__ = ('s', 'p')

    def __init__(self, raw):
        self.s = ''.join(format(b, '08b') for b in raw)
        self.p = 0

    def get(self, n):
        p = self.p
        self.p = p + n
        return int(self.s[p:p + n], 2) if n else 0

    def ints3(self, nbits, sizes):
        """receiveints: 3 integers packed in nbits bits as a mixed-radix number, least significant byte first."""
        big, shift = 0, 0
        while nbits > 8:
            big |= self.get(8) << shift
            shift += 8
            nbits -= 8
        if nbits > 0:
            big |= self.get(nbits) << shift
        c = big % sizes[2]
        big //= sizes[2]
        b = big % sizes[1]
        a = big // sizes[1]
        return a, b, c


def decode(data, off):
    natoms, step, time, box, lsize, p = _header(data, off)
    if lsize <= 9:
        v = struct.unpack_from(f'>{3 * lsize}f', data, p)
        return natoms, step, time, box, list(v[0::3]), list(v[1::3]), list(v[2::3])
    precision = struct.unpack_from('>f', data, p)[0]
    minint = struct.unpack_from('>3i', data, p + 4)
    maxint = struct.unpack_from('>3i', data, p + 16)
    smallidx = struct.unpack_from('>i', data, p + 28)[0]
    nbytes = struct.unpack_from('>i', data, p + 32)[0]
    bits = _Bits(data[p + 36:p + 36 + nbytes])

    sizeint = [maxint[k] - minint[k] + 1 for k in range(3)]
    if (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff:
        bitsizeint = [s.bit_length() for s in sizeint]
        bitsize = 0
    else:
        bitsizeint = None
        bitsize = (sizeint[0] * sizeint[1] * sizeint[2]).bit_length()
    smaller = MAGICINTS[max(FIRSTIDX, smallidx - 1)] // 2
    smallnum = MAGICINTS[smallidx] // 2
    sizesmall = [MAGICINTS[smallidx]] * 3
    inv = 1.0 / precision
    X, Y, Z = [], [], []
    i = 0
    run = 0
    while i < lsize:
        if bitsize == 0:
            t0, t1, t2 = bits.get(bitsizeint[0]), bits.get(bitsizeint[1]), bits.get(bitsizeint[2])
        else:
            t0, t1, t2 = bits.ints3(bitsize, sizeint)
        i += 1
        t0 += minint[0]
        t1 += minint[1]
        t2 += minint[2]
        p0, p1, p2 = t0, t1, t2
        flag = bits.get(1)
        is_smaller = 0
        if flag == 1:
            run = bits.get(5)
            is_smaller = run % 3
            run -= is_smaller
            is_smaller -= 1
        if run > 0:
            for k in range(0, run, 3):
                s0, s1, s2 = bits.ints3(smallidx, sizesmall)
                i += 1
                s0 += p0 - smallnum
                s1 += p1 - smallnum
                s2 += p2 - smallnum
                if k == 0:
                    # the first small atom is written before the large one (water trick)
                    s0, p0 = p0, s0
                    s1, p1 = p1, s1
                    s2, p2 = p2, s2
                    X.append(p0 * inv)
                    Y.append(p1 * inv)
                    Z.append(p2 * inv)
                else:
                    p0, p1, p2 = s0, s1, s2
                X.append(s0 * inv)
                Y.append(s1 * inv)
                Z.append(s2 * inv)
        else:
            X.append(t0 * inv)
            Y.append(t1 * inv)
            Z.append(t2 * inv)
        smallidx += is_smaller
        if is_smaller < 0:
            smallnum = smaller
            smaller = MAGICINTS[smallidx - 1] // 2 if smallidx > FIRSTIDX else 0
        elif is_smaller > 0:
            smaller = smallnum
            smallnum = MAGICINTS[smallidx] // 2
        sizesmall = [MAGICINTS[smallidx]] * 3
    if len(X) != lsize:
        raise ValueError(f'decoded {len(X)} atoms, expected {lsize}')
    return natoms, step, time, box, X, Y, Z
