
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


// integrators/mlt.cpp*
#include "integrators/pfmlt.h"
#include "integrators/bdpt.h"
#include "scene.h"
#include "film.h"
#include "integrator.h"
#include "camera.h"
#include "stats.h"
#include "filters/box.h"
#include "paramset.h"
#include "sampling.h"
#include "progressreporter.h"

namespace pbrt {

STAT_PERCENT("Integrator/Acceptance rate", acceptedpfMutations, totalpfMutations);

// MLTSampler Constants
static const int cameraStreamIndex = 0;
static const int lightStreamIndex = 1;
static const int connectionStreamIndex = 2;
static const int nSampleStreams = 3;

// PFMLT Method Definitions
Spectrum PFMLTIntegrator::L(const Scene &scene, MemoryArena &arena,
                            const std::unique_ptr<Distribution1D> &lightDistr,
                            const std::unordered_map<const Light *, size_t> &lightToIndex,
                            MLTSampler &sampler, int depth, Point2f *pRaster) {
    sampler.StartStream(cameraStreamIndex);
    // Determine the number of available strategies and pick a specific one
    int s, t, nStrategies;
    if (depth == 0) {
        nStrategies = 1;
        s = 0;
        t = 2;
    } else {
        nStrategies = depth + 2;
        s = std::min((int)(sampler.Get1D() * nStrategies), nStrategies - 1);
        t = nStrategies - s;
    }
    
    // Generate a camera subpath with exactly _t_ vertices
    Vertex *cameraVertices = arena.Alloc<Vertex>(t);
    Bounds2f sampleBounds = (Bounds2f)camera->film->GetSampleBounds();
    *pRaster = sampleBounds.Lerp(sampler.Get2D());
    if (GenerateCameraSubpath(scene, sampler, arena, t, *camera, *pRaster,
                              cameraVertices) != t)
        return Spectrum(0.f);
    
    // Generate a light subpath with exactly _s_ vertices
    sampler.StartStream(lightStreamIndex);
    Vertex *lightVertices = arena.Alloc<Vertex>(s);
    if (GenerateLightSubpath(scene, sampler, arena, s, cameraVertices[0].time(),
                             *lightDistr, lightToIndex, lightVertices) != s)
        return Spectrum(0.f);
    
    // Execute connection strategy and return the radiance estimate
    sampler.StartStream(connectionStreamIndex);
    return ConnectBDPT(scene, lightVertices, cameraVertices, s, t, *lightDistr,
                       lightToIndex, *camera, sampler, pRaster) *
    nStrategies;
}

void PFMLTIntegrator::Render(const Scene &scene) {
    std::unique_ptr<Distribution1D> lightDistr =
    ComputeLightPowerDistribution(scene);
    
    // Compute a reverse mapping from light pointers to offsets into the
    // scene lights vector (and, equivalently, offsets into
    // lightDistr). Added after book text was finalized; this is critical
    // to reasonable performance with 100s+ of light sources.
    std::unordered_map<const Light *, size_t> lightToIndex;
    for (size_t i = 0; i < scene.lights.size(); ++i)
        lightToIndex[scene.lights[i].get()] = i;
    
    // Generate bootstrap samples and compute normalization constant $b$
    int nBootstrapSamples = nBootstrap * (maxDepth + 1);
    std::vector<Float> bootstrapWeights(nBootstrapSamples, 0);
    if (scene.lights.size() > 0) {
        ProgressReporter progress(nBootstrap / 256,
                                  "Generating bootstrap paths");
        std::vector<MemoryArena> bootstrapThreadArenas(MaxThreadIndex());
        int chunkSize = Clamp(nBootstrap / 128, 1, 8192);
        ParallelFor([&](int i) {
            // Generate _i_th bootstrap sample
            MemoryArena &arena = bootstrapThreadArenas[ThreadIndex];
            for (int depth = 0; depth <= maxDepth; ++depth) {
                int rngIndex = i * (maxDepth + 1) + depth;
                MLTSampler sampler(mutationsPerPixel, rngIndex, sigma,
                                   largeStepProbability, nSampleStreams);
                Point2f pRaster;
                bootstrapWeights[rngIndex] =
                L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pRaster).y();
                arena.Reset();
            }
            if ((i + 1 % 256) == 0) progress.Update();
        }, nBootstrap, chunkSize);
        progress.Done();
    }
    Distribution1D bootstrap(&bootstrapWeights[0], nBootstrapSamples);
    Float b = bootstrap.funcInt * (maxDepth + 1);
    
#define ZLB3 1
    int totalM=0, acceptM=0, againM=0, zeroM=0;//zlb3
    
    // Run _nChains_ Markov chains in parallel
    Film &film = *camera->film;
    int64_t nTotalMutations =
    (int64_t)mutationsPerPixel * (int64_t)film.GetSampleBounds().Area();
    if (scene.lights.size() > 0) {
        const int progressFrequency = 32768;
        ProgressReporter progress(nTotalMutations / progressFrequency,
                                  "Rendering");
        
        ParallelFor([&](int i) {
            int64_t nChainMutations =
            std::min((i + 1) * nTotalMutations / nChains, nTotalMutations) -
            i * nTotalMutations / nChains;
            // Follow {i}th Markov chain for _nChainMutations_
            MemoryArena arena;
            
            // Select initial state from the set of bootstrap samples
            RNG rng(i);
            int bootstrapIndex;
            bootstrapIndex= bootstrap.SampleDiscrete(rng.UniformFloat());
            
            
            int depth = bootstrapIndex % (maxDepth + 1);
            
            // Initialize local variables for selected state
            MLTSampler sampler(mutationsPerPixel, bootstrapIndex, sigma,
                               largeStepProbability, nSampleStreams);
            Point2f pCurrent;
            Spectrum LCurrent =
            L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pCurrent);//zlb2
#if ZLB3
            int CNT = 0;
            int MAXREP;
            if(maxRep == 0){
                MAXREP = depth + 2;
            }
            else{
                MAXREP = std::min(maxRep, depth + 2);
            }
#endif//ZLB3
            // Run the Markov chain for _nChainMutations_ steps
            for (int64_t j = 0; j < nChainMutations; ++j) {
                sampler.StartIteration();
                Point2f pProposed;
                Spectrum LProposed;
                LProposed = L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pProposed);
                
                // Compute acceptance probability for proposed sample
                Float accept = std::min((Float)1, LProposed.y() / LCurrent.y());
                
                totalM++;//zlb3
#if ZLB3//just try again
                if(accept < 0.000001 && CNT < MAXREP){
                    CNT ++;
                    sampler.Reject();
//                    j --;
                    againM++;//zlb3
                    continue;
                }
                CNT = 0;
#endif//ZLB3
                if(accept < 0.000001){
                    zeroM++;
                }
                
                // Splat both current and proposed samples to _film_
                if (accept > 0)
                    film.AddSplat(pProposed,
                                  LProposed * accept / LProposed.y());
                film.AddSplat(pCurrent, LCurrent * (1 - accept) / LCurrent.y());
                
                // Accept or reject the proposal
                if (rng.UniformFloat() < accept) {
                    pCurrent = pProposed;
                    LCurrent = LProposed;
                    acceptM++;//zlb3
                    sampler.Accept();
                    ++acceptedpfMutations;
                } else
                    sampler.Reject();
                ++totalpfMutations;
                if ((i * nTotalMutations / nChains + j) % progressFrequency ==
                    0)
                    progress.Update();
                arena.Reset();
            }
        }, nChains);
        progress.Done();
    }
    
    printf("totalM=%d, acceptM=%d, againM=%d, zeroM=%d\n", totalM, acceptM, againM, zeroM);//zlb3
    
    // Store final image computed with MLT
    camera->film->WriteImage(b / mutationsPerPixel);
}

PFMLTIntegrator *CreatePFMLTIntegrator(const ParamSet &params,
                                       std::shared_ptr<const Camera> camera) {
    int maxRep = params.FindOneInt("maxrep", 0);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
    int64_t nChains = params.FindOneInt("chains", 1000);
    int mutationsPerPixel = params.FindOneInt("mutationsperpixel", 100);
    Float largeStepProbability =
    params.FindOneFloat("largestepprobability", 0.3f);
    Float sigma = params.FindOneFloat("sigma", .01f);
    if (PbrtOptions.quickRender) {
        mutationsPerPixel = std::max(1, mutationsPerPixel / 16);
        nBootstrap = std::max(1, nBootstrap / 16);
    }
    return new PFMLTIntegrator(camera, maxRep, maxDepth, nBootstrap, nChains,
                               mutationsPerPixel, sigma, largeStepProbability);
}

}  // namespace pbrt
