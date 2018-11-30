
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// media/grid.cpp*
#include "media/gridrm.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

namespace pbrt {

STAT_RATIO("Media/Grid steps per Tr() call", nTrSteps, nTrCalls);

// GridDensityMedium Method Definitions
Float GridDensityMediumRayMarch::Density(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // Trilinearly interpolate density values to compute local density
    Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
    Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
    Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
    Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
    Float d0 = Lerp(d.y, d00, d10);
    Float d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

Spectrum GridDensityMediumRayMarch::Sample(const Ray &rWorld, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // Run delta-tracking iterations to sample a medium interaction
    Float t = tMin;

    // sample an optical depth 
    Float sampledTau = -std::log(1 - sampler.Get1D());
    Float accumTau = 0.f;
    // ray march till we accrue optical depth
    while (true) {
        if (t >= tMax) break;
        // Density returns values [0-1]
        
        Float temp = sigma_t * Density(ray(t)) * stepSize;
        Float tau = accumTau + temp;
        if (tau > sampledTau) {
            // Populate _mi_ with medium interaction information and return
            PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
            *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                                    phase);
            return sigma_s / sigma_t;
        }
        else
        {
            // proceed further
            accumTau += temp;
            t += stepSize;
        }
    }
    return Spectrum(1.f);
}

Spectrum GridDensityMediumRayMarch::Tr(const Ray &rWorld, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    ++nTrCalls;

    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);
    // Perform ray marching to get the optical depth and then exponentiate
    Float t = tMin;
    Float tau = 0.f;
    while (true) {
        ++nTrSteps;
        t += stepSize;
        if (t >= tMax) break;
        tau += sigma_t * Density(ray(t)) * stepSize;
    }

    // compute transmittance
    Float Tr = std::expf(-tau);
    return Spectrum(Tr);
}

}  // namespace pbrt
