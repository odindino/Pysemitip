# Detailed Mapping: potexpand-6.1.f ↔ physics/core/potential.py

## Overview
**Fortran File**: `src/fortran/MultInt/potexpand-6.1.f`  
**Python Module**: `src/physics/core/potential.py::PotentialProcessor.expand_potential()`  
**Completion Status**: 45% - Partial implementation, missing image potential and detailed expansion algorithms  
**Priority**: MEDIUM - Important for accurate Schrödinger equation integration

## Function Purpose
POTEXPAND expands the grid of z-values and corresponding potentials in vacuum and semiconductor regions to create approximately equal z-spacing suitable for integrating the Schrödinger equation. This is crucial for accurate tunneling current calculations.

## Fortran Code Structure
```fortran
SUBROUTINE POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NS,NSDIM,BARR,NBARR1,
    BARR2,NBARR2,NVDIM1,NVDIM2,PROF,PROF2,NSDIM2,S2,NS2,VACSTEP,
    SEMSTEP,JSEM,NEXSEM,NEXVAC,IWRIT)
```

## Line-by-Line Mapping

### Input Parameters (Lines 12-38)
| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 12-15: IMPOT | Not implemented | ❌ | Image potential flag missing |
| 15: SEP | method parameters | ⚠️ | Tip-sample separation |
| 16: NV | grid.params.nv | ✅ | Vacuum grid points |
| 17: Pot0p | method parameters | ⚠️ | Surface potential |
| 18-19: S array | grid.zs | ✅ | Semiconductor z-coordinates |
| 20-23: BARR arrays | method parameters | ⚠️ | Vacuum potential arrays |
| 24-27: PROF arrays | method parameters | ⚠️ | Semiconductor potential arrays |
| 28-33: Grid parameters | expansion_factor | ⚠️ | Expansion control parameters |

### Vacuum Expansion Section (Lines 48-84)

#### Expansion Factor Calculation (Lines 48-56)
```fortran
nexpan=MAX0(1,NINT((SEP/NV)/VACSTEP))
IF (IMPOT.EQ.1) NEXPAN=NEXPAN*10
NEXVAC=NEXPAN
```

**Python Implementation:**
```python
# In expand_potential() - missing detailed vacuum expansion
delv_new = old_params.delv / expansion_factor
nv_new = (old_params.nv - 1) * expansion_factor + 1
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 48-49 | delv_new calculation | ⚠️ | Simple division vs adaptive expansion |
| 52-53 | Not implemented | ❌ | Image potential 10x expansion missing |
| 54 | expansion_factor | ✅ | Store expansion factor |

#### Linear Interpolation (Lines 57-67)
```fortran
NBARR2=nexpan*(NBARR1-1)+1
BARR2(NBARR2)=BARR(NBARR1)
DO 150 J=NBARR1-1,1,-1
   B2=BARR(J+1)
   B1=BARR(J)
   DO 140 K=nexpan-1,0,-1
      BARR2((J-1)*nexpan+K+1)=(B2*FLOAT(K)+B1*FLOAT(nexpan-K))/nexpan
```

**Python Implementation:**
```python
# In expand_potential() - uses scipy interpolation
interp_func = interpolate.RegularGridInterpolator(
    (r_old, z_old, phi_old), potential_3d, method='linear')
potential_new = interp_func(points).reshape(nr_new, len(z_new), np_new)
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 57 | nv_new calculation | ✅ | New array size |
| 58-67 | scipy interpolation | ⚠️ | Different method but equivalent |

#### Image Potential Correction (Lines 74-80)
```fortran
lambda=3.81**2*0.1*alog(2.)/(2.*2.*sep)
IF (IMPOT.EQ.1) THEN
do 200 j=2,NBARR2-1
   barr2(j)=barr2(j)-1.15*lambda*(NBARR2-1.)**2/
            ((j-1.)*(float(NBARR2)-j))
```

**Python Implementation:**
```python
# NOT IMPLEMENTED - critical missing feature
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 74 | Not implemented | ❌ | Lambda calculation missing |
| 75-80 | Not implemented | ❌ | Image potential correction missing |

### Semiconductor Expansion Section (Lines 89-144)

#### Initialize Expansion Array (Lines 89-95)
```fortran
DO 300 J=1,NS
   NEXSEM(J)=0
300 CONTINUE
nexpan=max0(1,NINT(2.*S(1)/SEMSTEP))
```

**Python Implementation:**
```python
# In expand_potential() - simplified
dels_new = old_params.dels / expansion_factor
ns_new = (old_params.ns - 1) * expansion_factor + 1
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 89-91 | Not needed | ✅ | Python handles automatically |
| 94 | expansion_factor logic | ⚠️ | Different calculation method |

#### Adaptive Expansion Loop (Lines 101-144)
```fortran
DO 570 J=1,NS
   IF (J.EQ.1) THEN
      NEXPAN=MAX0(1,NINT(S(1)/SEMSTEP))
   ELSE
      NEXPAN=MAX0(1,NINT((S(J)-S(J-1))/SEMSTEP))
   END IF
   IF (MOD(NEXPAN,2).EQ.0) NEXPAN=NEXPAN+1
```

**Python Implementation:**
```python
# Uses uniform interpolation instead of adaptive
z_new = np.concatenate([
    -np.linspace(0, old_params.smax, ns_new)[::-1][:-1],
    np.linspace(0, old_params.vmax, nv_new)
])
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 101-108 | Uniform spacing | ⚠️ | Python uses uniform vs adaptive |
| 109 | Not implemented | ❌ | Odd number enforcement missing |

#### Point Assignment and Interpolation (Lines 110-138)
```fortran
DO 560 K=1,NEXPAN
   KK=KK+1
   IF (J.EQ.1) THEN
      JSEM(KK)=J
   ELSE
      IF (K.LE.(NEXPAN/2)) THEN
         JSEM(KK)=J-1
      ELSE
         JSEM(KK)=J
      END IF
   END IF
   IF (J.EQ.1) THEN
      PROF2(KK)=((NEXPAN-K)*Pot0p+(K)*PROF(J))/FLOAT(NEXPAN)
   ELSE
      PROF2(KK)=((NEXPAN-K)*PROF(J-1)+(K)*PROF(J))/FLOAT(NEXPAN)
   END IF
```

**Python Implementation:**
```python
# Scipy handles interpolation automatically
potential_new = interp_func(points).reshape(nr_new, len(z_new), np_new)
```

| Fortran Lines | Python Equivalent | Status | Notes |
|---------------|-------------------|---------|-------|
| 110-120 | Automatic indexing | ✅ | Python handles automatically |
| 126-132 | scipy interpolation | ✅ | Different method but equivalent |
| 134-137 | scipy interpolation | ✅ | Z-coordinate interpolation |

## Python Implementation Analysis

### Current Implementation
```python
def expand_potential(self, potential_3d: np.ndarray,
                   expansion_factor: int = 2) -> Tuple[np.ndarray, Grid3D]:
    """Expand potential to finer grid."""
    # Create new grid parameters
    old_params = self.grid.params
    
    # New grid dimensions
    nr_new = (old_params.nr - 1) * expansion_factor + 1
    nv_new = (old_params.nv - 1) * expansion_factor + 1
    ns_new = (old_params.ns - 1) * expansion_factor + 1
    
    # Use scipy interpolation
    interp_func = interpolate.RegularGridInterpolator(
        (r_old, z_old, phi_old), potential_3d, method='linear')
    potential_new = interp_func(points).reshape(nr_new, len(z_new), np_new)
    
    return potential_new, None
```

### Fortran Algorithm (Simplified)
```fortran
! 1. Calculate vacuum expansion factor
nexpan_vac = MAX(1, NINT((SEP/NV)/VACSTEP))
if (IMPOT == 1) nexpan_vac = nexpan_vac * 10

! 2. Expand vacuum potential with linear interpolation
DO J=NBARR1-1,1,-1
   DO K=nexpan-1,0,-1
      BARR2(index) = (B2*K + B1*(nexpan-K))/nexpan

! 3. Apply image potential correction
if (IMPOT == 1) then
   lambda = 3.81^2 * 0.1 * log(2) / (4 * sep)
   BARR2(j) = BARR2(j) - 1.15*lambda*(NBARR2-1)^2/((j-1)*(NBARR2-j))

! 4. Expand semiconductor with adaptive spacing
DO J=1,NS
   nexpan_local = MAX(1, NINT(delta_z/SEMSTEP))
   if (MOD(nexpan_local,2) == 0) nexpan_local = nexpan_local + 1
```

## Key Differences

### 1. Expansion Strategy
- **Fortran**: Adaptive expansion based on local grid spacing
- **Python**: Uniform expansion with fixed factor

### 2. Image Potential
- **Fortran**: Full image potential correction with lambda calculation
- **Python**: Not implemented

### 3. Interpolation Method
- **Fortran**: Custom linear interpolation with exact point assignment
- **Python**: Scipy RegularGridInterpolator

### 4. Grid Spacing Control
- **Fortran**: Enforces odd number of expansion points, adaptive local spacing
- **Python**: Uniform spacing throughout

## Missing Features in Python

### Critical Missing Features
1. **Image Potential Correction** (Lines 74-80)
   - Lambda calculation: `3.81²×0.1×ln(2)/(4×sep)`
   - Correction formula: `-1.15×λ×(N-1)²/((j-1)×(N-j))`

2. **Adaptive Spacing Algorithm** (Lines 101-108)
   - Local expansion factor based on grid spacing
   - Odd number enforcement for numerical stability

3. **Point Assignment Logic** (Lines 110-120)
   - JSEM array for tracking correspondence between grids
   - Half-interval assignment for accurate interpolation

### Implementation Gaps
1. **Surface potential handling** (Pot0p parameter)
2. **Expansion factor arrays** (NEXSEM, NEXVAC)
3. **Debug output control** (IWRIT parameter)

## Usage Context

### In Fortran INTCURR
```fortran
CALL POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NSP,NSDIM,BARR,NBARR1,BARR2,
    NBARR2,NVDIM1,NVDIM2,CBPROF,PROF2,NSDIM2,S2,NS2,VACSTEP,SEMSTEP,
    JSEM,NEXSEM,NEXVAC,IWRIT)
```

### In Python SchrodingerSolver
```python
# Currently not used - integration needed
expanded_potential, expanded_grid = processor.expand_potential(
    potential_3d, expansion_factor=4)
```

## Integration Status

### Current Integration: 45%
- ✅ Basic grid expansion framework
- ✅ Scipy interpolation infrastructure  
- ⚠️ Simple uniform expansion (not adaptive)
- ❌ Image potential correction missing
- ❌ Adaptive spacing algorithm missing
- ❌ Point correspondence tracking missing

### Required for 95% Completion
1. Implement image potential correction
2. Add adaptive expansion algorithm
3. Include point assignment tracking
4. Add surface potential handling
5. Integrate with SchrodingerSolver

## Priority Assessment
**MEDIUM Priority** - Important for accurate current calculations but not critical for basic functionality. The current scipy-based implementation provides reasonable interpolation, but lacks the physical corrections (image potential) and adaptive spacing that make the Fortran version more accurate for tunneling calculations.
