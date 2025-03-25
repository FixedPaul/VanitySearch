/*
 * This file is part of the VanitySearch distribution (https://github.com/JeanLucPons/VanitySearch).
 * Copyright (c) 2019 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// CUDA Kernel main function
// Compute SecpK1 keys and calculate RIPEMD160(SHA256(key)) then check prefix
// For the kernel, we use a 16 bits prefix lookup table which correspond to ~3 Base58 characters
// A second level lookup table contains 32 bits prefix (if used)
// (The CPU computes the full address and check the full prefix)
//
// We use affine coordinates for elliptic curve point (ie Z=1)

__device__ __noinline__ void CheckPointPattern(uint32_t *_h, int32_t incr, int32_t endo, int32_t mode,prefix_t *prefix,
                                        uint32_t *lookup32, uint32_t maxFound, uint32_t *out,int type) {

  uint32_t   pos;  
  char       add[48];

  // No lookup compute address and return
  char* pattern = (char*)lookup32;
  _GetAddress(type, _h, add);
  if (_Match(add, pattern)) {
      // found
      pos = atomicAdd(out, 1);
      if (pos < maxFound) {
          uint32_t   tid = (blockIdx.x * blockDim.x) + threadIdx.x;
          out[pos * ITEM_SIZE32 + 1] = tid;
          out[pos * ITEM_SIZE32 + 2] = (uint32_t)(incr << 16) | (uint32_t)(mode << 15) | (uint32_t)(endo);
          out[pos * ITEM_SIZE32 + 3] = _h[0];
          out[pos * ITEM_SIZE32 + 4] = _h[1];
          out[pos * ITEM_SIZE32 + 5] = _h[2];
          out[pos * ITEM_SIZE32 + 6] = _h[3];
          out[pos * ITEM_SIZE32 + 7] = _h[4];
      }
  }

}


__device__ __noinline__ void CheckPointLUT(uint32_t* _h, int32_t incr, int32_t endo, int32_t mode, prefix_t* prefix,
    uint32_t* lookup32, uint32_t maxFound, uint32_t* out, int type) {

    uint32_t   off;
    prefixl_t  l32;
    prefix_t   pr0;
    prefix_t   hit;
    uint32_t   pos;
    uint32_t   st;
    uint32_t   ed;
    uint32_t   mi;
    uint32_t   lmi;

    // Lookup table
    pr0 = *(prefix_t*)(_h);
    hit = prefix[pr0];

    if (hit) {

        if (lookup32) {
            off = lookup32[pr0];
            l32 = _h[0];
            st = off;
            ed = off + hit - 1;
            while (st <= ed) {
                mi = (st + ed) / 2;
                lmi = lookup32[mi];
                if (l32 < lmi) {
                    ed = mi - 1;
                }
                else if (l32 == lmi) {
                    // found
                    goto addItem;
                }
                else {
                    st = mi + 1;
                }
            }
            return;
        }

    addItem:

        pos = atomicAdd(out, 1);
        if (pos < maxFound) {
            uint32_t   tid = (blockIdx.x * blockDim.x) + threadIdx.x;
            out[pos * ITEM_SIZE32 + 1] = tid;
            out[pos * ITEM_SIZE32 + 2] = (uint32_t)(incr << 16) | (uint32_t)(mode << 15) | (uint32_t)(endo);
            out[pos * ITEM_SIZE32 + 3] = _h[0];
            out[pos * ITEM_SIZE32 + 4] = _h[1];
            out[pos * ITEM_SIZE32 + 5] = _h[2];
            out[pos * ITEM_SIZE32 + 6] = _h[3];
            out[pos * ITEM_SIZE32 + 7] = _h[4];
        }

    }
        


}

// -----------------------------------------------------------------------------------------  P2PKH
 

#define CHECK_P2PKH_COMP(_incr) {                                             \
_GetHash160Comp(px, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(px, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                        \
CheckPointLUT(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                        \
CheckPointLUT(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                        \
CheckPointLUT(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysComp(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    CHECK_P2PKH_COMP(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);


        CHECK_P2PKH_COMP(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        CHECK_P2PKH_COMP(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    CHECK_P2PKH_COMP(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

// -----------------------------------------------------------------------------------------

#define CHECK_P2PKH_UNCOMP(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160UnComp(px, py, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointLUT(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysUnComp(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {


    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2PKH_UNCOMP(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2PKH_UNCOMP(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2PKH_UNCOMP(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2PKH_UNCOMP(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

#define CHECK_P2PKH_BOTH(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160UnComp(px, py, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(px, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(px, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointLUT(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysBoth(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2PKH_BOTH(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2PKH_BOTH(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2PKH_BOTH(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2PKH_BOTH(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}


// -----------------------------------------------------------------------------------------  P2SH

#define CHECK_P2SH_COMP(_incr) {                                             \
_GetHash160P2SHComp(px, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(px, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                        \
CheckPointLUT(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                        \
CheckPointLUT(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                        \
CheckPointLUT(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysCompP2SH(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    CHECK_P2SH_COMP(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);


        CHECK_P2SH_COMP(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        CHECK_P2SH_COMP(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    CHECK_P2SH_COMP(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

// -----------------------------------------------------------------------------------------

#define CHECK_P2SH_UNCOMP(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160P2SHUnComp(px, py, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointLUT(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysUnCompP2SH(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {


    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2SH_UNCOMP(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2SH_UNCOMP(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2SH_UNCOMP(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2SH_UNCOMP(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

#define CHECK_P2SH_BOTH(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160P2SHUnComp(px, py, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(px, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(px, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointLUT(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointLUT(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointLUT(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                           \
CheckPointLUT(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                            \
CheckPointLUT(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysBothP2SH(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2SH_BOTH(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2SH_BOTH(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2SH_BOTH(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2SH_BOTH(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}


//////////////////////////////////////////////////////////////////////////////////// PATTERN SEARCH


// -----------------------------------------------------------------------------------------  P2PKH


#define CHECK_P2PKH_COMP_PATTERN(_incr) {                                             \
_GetHash160Comp(px, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(px, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                        \
CheckPointPattern(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                        \
CheckPointPattern(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                        \
CheckPointPattern(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysCompPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    CHECK_P2PKH_COMP_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);


        CHECK_P2PKH_COMP_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        CHECK_P2PKH_COMP_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    CHECK_P2PKH_COMP_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

// -----------------------------------------------------------------------------------------

#define CHECK_P2PKH_UNCOMP_PATTERN(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160UnComp(px, py, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointPattern(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysUnCompPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {


    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2PKH_UNCOMP_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2PKH_UNCOMP_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2PKH_UNCOMP_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2PKH_UNCOMP_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

#define CHECK_P2PKH_BOTH_PATTERN(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160UnComp(px, py, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(px, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(px, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160UnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160UnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointPattern(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2PKH);    \
_GetHash160Comp(pex, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);     \
_GetHash160Comp(pex, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2PKH);    \
}

__device__ void ComputeKeysBothPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2PKH_BOTH_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2PKH_BOTH_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2PKH_BOTH_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2PKH_BOTH_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}


// -----------------------------------------------------------------------------------------  P2SH

#define CHECK_P2SH_COMP_PATTERN(_incr) {                                             \
_GetHash160P2SHComp(px, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(px, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                        \
CheckPointPattern(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                        \
CheckPointPattern(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                        \
CheckPointPattern(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysCompP2SHPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    CHECK_P2SH_COMP_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);


        CHECK_P2SH_COMP_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        CHECK_P2SH_COMP_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    CHECK_P2SH_COMP_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

// -----------------------------------------------------------------------------------------

#define CHECK_P2SH_UNCOMP_PATTERN(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160P2SHUnComp(px, py, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointPattern(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysUnCompP2SHPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {


    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2SH_UNCOMP_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2SH_UNCOMP_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2SH_UNCOMP_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2SH_UNCOMP_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}

#define CHECK_P2SH_BOTH_PATTERN(_incr) {                                             \
ModNeg256(pyn, py);                                                             \
_GetHash160P2SHUnComp(px, py, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(px, pyn, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(px, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(px, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 0, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta);                                                     \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t*)h);                                  \
CheckPointPattern(h, -(_incr), 1, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 1, true, sPrefix, lookup32, maxFound, out, P2SH);    \
_ModMult(pex, px, _beta2);                                                    \
_GetHash160P2SHUnComp(pex, py, (uint8_t *)h);                        \
CheckPointPattern(h, (_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHUnComp(pex, pyn, (uint8_t *)h);                        \
CheckPointPattern(h, -(_incr), 2, false, sPrefix, lookup32, maxFound, out, P2SH);    \
_GetHash160P2SHComp(pex, 0, (uint8_t *)h);                                           \
CheckPointPattern(h, (_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);     \
_GetHash160P2SHComp(pex, 1, (uint8_t *)h);                                            \
CheckPointPattern(h, -(_incr), 2, true, sPrefix, lookup32, maxFound, out, P2SH);    \
}

__device__ void ComputeKeysBothP2SHPattern(uint64_t* startx, uint64_t* starty, prefix_t* sPrefix, uint32_t* lookup32, uint32_t maxFound, uint32_t* out) {

    uint64_t dx[4];
    uint64_t px[4];
    uint64_t py[4];
    uint64_t dy[4];
    uint64_t sxn[4];
    uint64_t syn[4];
    uint64_t sx[4];
    uint64_t sy[4];
    uint64_t sx_gx[4];
    uint64_t inverse[5];

    uint32_t   h[5];
    uint64_t   pex[4];
    uint64_t   pyn[4];

    uint64_t subp[GRP_SIZE / 2][4];


    __syncthreads();
    Load256A(sx, startx);
    Load256A(sy, starty);


    uint32_t i;

    // Check starting point
    Load256(px, sx);
    Load256(py, sy);
    CHECK_P2SH_BOTH_PATTERN(GRP_SIZE / 2);

    __syncthreads();

    ModSub256(sxn, _2Gnx, sx);
    Load256(subp[GRP_SIZE / 2 - 1], sxn);
    for (i = GRP_SIZE / 2 - 1; i > 0; i--) {
        ModSub256(syn, Gx[i], sx);
        _ModMult(sxn, syn);
        Load256(subp[i - 1], sxn);
    }

    ModSub256(inverse, Gx[0], sx);
    _ModMult(inverse, sxn);


    inverse[4] = 0;
    _ModInv(inverse);

    __syncthreads();

    ModNeg256(syn, sy);
    ModNeg256(sxn, sx);

    for (i = 0; i < GRP_SIZE / 2 - 1; i++) {

        __syncthreads();
        ModSub256(sx_gx, Gx[i], sxn);

        _ModMult(dx, subp[i], inverse);

        //////////////////

        ModSub256(dy, Gy[i], sy);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, sx, px);
        _ModMult(py, dy);
        ModSub256(py, sy);

        CHECK_P2SH_BOTH_PATTERN(GRP_SIZE / 2 + (i + 1));

        //////////////////

        __syncthreads();

        ModSub256(dy, syn, Gy[i]);
        _ModMult(dy, dx);
        _ModSqr(px, dy);
        ModSub256(px, sx_gx);

        ModSub256(py, px, sx);
        _ModMult(py, dy);
        ModSub256(py, syn, py);

        CHECK_P2SH_BOTH_PATTERN(GRP_SIZE / 2 - (i + 1));

        //////////////////

        ModSub256(dx, Gx[i], sx);
        _ModMult(inverse, dx);

    }

    __syncthreads();

    _ModMult(dx, subp[i], inverse);

    ModSub256(dy, syn, Gy[i]);
    _ModMult(dy, dx);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, Gx[i]);

    ModSub256(py, px, sx);
    _ModMult(py, dy);
    ModSub256(py, syn, py);

    CHECK_P2SH_BOTH_PATTERN(0);

    //////////////////

    __syncthreads();

    ModSub256(dy, _2Gny, sy);
    ModSub256(dx, Gx[i], sx);
    _ModMult(inverse, dx);

    _ModMult(dy, inverse);
    _ModSqr(px, dy);
    ModSub256(px, sx);
    ModSub256(px, _2Gnx);

    ModSub256(py, _2Gnx, px);
    _ModMult(py, dy);
    ModSub256(py, _2Gny);


    __syncthreads();
    Store256A(startx, px);
    Store256A(starty, py);

}