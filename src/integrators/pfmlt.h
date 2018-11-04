
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PFMLT_H
#define PBRT_INTEGRATORS_PFMLT_H

// integrators/mlt.h*
#include "pbrt.h"
#include "integrator.h"
#include "mlt.h"
#include "spectrum.h"
#include "film.h"
#include "rng.h"
#include <unordered_map>

namespace pbrt {
    
    // MLT Declarations
    class PFMLTIntegrator : public Integrator {
    public:
        // MLTIntegrator Public Methods
        PFMLTIntegrator(std::shared_ptr<const Camera> camera, int maxRep, int maxDepth,
                        int nBootstrap, int nChains, int mutationsPerPixel,
                        Float sigma, Float largeStepProbability)
        : camera(camera),
        maxRep(maxRep),
        maxDepth(maxDepth),
        nBootstrap(nBootstrap),
        nChains(nChains),
        mutationsPerPixel(mutationsPerPixel),
        sigma(sigma),
        largeStepProbability(largeStepProbability) {}
        void Render(const Scene &scene);
        Spectrum L(const Scene &scene, MemoryArena &arena,
                   const std::unique_ptr<Distribution1D> &lightDistr,
                   const std::unordered_map<const Light *, size_t> &lightToIndex,
                   MLTSampler &sampler, int k, Point2f *pRaster);
        
    private:
        // MLTIntegrator Private Data
        std::shared_ptr<const Camera> camera;
        const int maxRep;
        const int maxDepth;
        const int nBootstrap;
        const int nChains;
        const int mutationsPerPixel;
        const Float sigma, largeStepProbability;
    };
    
    PFMLTIntegrator *CreatePFMLTIntegrator(const ParamSet &params,
                                           std::shared_ptr<const Camera> camera);
    
}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_PFMLT_H
