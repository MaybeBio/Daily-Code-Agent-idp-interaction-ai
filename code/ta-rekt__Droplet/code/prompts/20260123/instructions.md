# Overview

This is a design and iteration session to put together an agentic pipeline for experimental design, execution and data analysis of numerical experiments using this repository.

We are going to do some design iterations to describe systematically what this pipeline should look like, one section at a time. The objective for this is to be able to parallelize setting up experiments and working on new features to significantly organize, accelerate and improve the science which this code supports.

At any moment during this iteration phase, or during the actual design of an experiment, you can ask me questions about something you don't understand or raise and issue of you think something is confusing and can be clarified.

Let's start with some general context about the project.

# Scientific background on the Droplet project

This repository is part of a PhD project on the sensing capabilities of biomolecular condensates in the context of gene regulation. Its purpose is to simulate polymers interacting with bulk proteins with a lattice-gas model. The scientific goal is to understand whether condensate formation and polymer dynamics can make a promoter sense the concentration of a particular protein more sensitively than by diffusion alone. By sensing we mean how the statistics of the transcription process change as a function of concentration, and how informative these changes can be to a downstream observer of transcription process. We are particularly paying attention to the limits of sensing (in the sense of Fisher information rate) imposed by physics on the binding and transcription process of a promoter.


# Description of the existing pipeline

## Specification file (located in `spec`)

Experiments typically consist of Monte-Carlo simulations of polymers interacting with bulk species in various scenarios, varying some parameter and monitoring some output. Each experiment is described in a specification file, located in the spec directory. Here is the general structure of such files:

### physical elements of the simulation

Instantiates the lattice and each species to feature in the simulation (no initialization or parameter assignment). 

### simulation parameters

Variables relating to the simulation itself (total runtime, burn_in time, number of steps between each recorded frame). Which parameter is swept over, its values, number of replicas (repeats)

### target function

Sets up the function to be called at each point in parameter space. This function assigns to the different species the parameters they should take at this point in the sweep. The main loop for moving each species is contained here, and the desired data to be collected is described after this loop. A dictionary called res is filled up and returned.

## Sweep script (located in `sim/sweep.jl`)

This script is meant to be executed in parallel over a number of CPUs which is set upstream, in the `runJL` bash script which is described next. It takes in 4 command line arguments: 1. the name of the specification file, 2. the key of the cpu on which to run it, 3. the total number of cpus, and 4. the output directory to which results are written. It makes an `include()` to load the specification file passed from the command line, then runs a loop over the specified set of parameters, builds up the results into a dict and writes it to the output path.

## runJL bash script (located in the root directory)

This is the high level script which orchestrates parallel cpu execution of the parameter sweeps with a unix pipe.

## Output data organization

Data resulting from this is stored in the `data` directory, inside a folder named as the date of the experiment in the YYYYMMDD format, in a directory named after the specification file associated with the experiment, with a short indication of the specific parameter range, name, or values used in this particular run, just to distinguish it from previous experiments.

## Example usage

This bash command runs the `sim/sweep.jl` script, with the `prewet_transcription.jl` specification file on 16 CPUs in parallel, and stores the output in `data/20251218/prewet_transcription_32_1300-4000`:

`./runJL sim/sweep.jl spec/prewet_transcription.jl 16 data/20251218/prewet_transcription_32_1300-4000`


# Proposal for an agentic pipeline

I would like to build on this existing workflow to make it possible for AI agents to handle it autonomously. I would like us to start iterating on each of the following workflows separately, starting with experiment design.


## Desirata of the pipeline to be designed

First, we are going to develop an experimental design workflow where I can describe in detail the experiment I need to set up, and have you, the agent, make a detailed proposal, on which we can iterate together until I approve it. the experiment proposal format needs to be specified, it must consist at mimimum of: creating a spec file, a summary of the experiment, and a description of the resulting plots, expected results, the technical specs such as number of cores to run on, expected simulation time, varied parameters, type of output data, etc. After the proposal is validated, possibly after a few iterations, the experiment moves to execution.

Second, we will design the experimental execution pipeline, which takes the proposal and specification files, runs the required simulations, retrieves and organizes the output data, and runs validation tests to ensure that everything went smoothly.

Third, the data summarization and visualization pipeline will take output data, create a detailed report and produce the plots described in the proposal, or, if the experiment failed, come up with an explanation or hypothesis of what went wrong and propose steps to address it in a new iteration. If all this goes well, the outcome should be a polished and presentable report for this particular experiment.

### General requirements:

- Make as much use of traditional software tools as possible in order to minimize premium requests to github copilot. If anything can be done with a tool and without agent intervention, we should build it. We want to keep agent intervention out of this workflow for all tasks that don't require intelligent analysis, hypothesis generation and testing, summarization and debugging (e.g. anything that can be automated with scripts).

- Create scripts for user notifications at various points in the lifecycle of an experiment: when a simulation starts and finishes executing, when data has been pulled into local folders, final report ready, signaling a problem, etc.

- Depart as little as possible from the existing structure and file content of the repo. Add as little context and harnesses as required for agents to navigate and work effectively, keeping the repo as it is to the maximum possible extent. It should remain a human-readable repo, because I am the main developer and I should understand each change being introduced.

# Objective

I want you to take a look at the repo, understand its structure, ask me questions if there's something you don't understand. Then, start by telling me what you think is the best way to set up context and scaffolding to make agents turn detailed descriptions of experiments into specification files. We will brainstorm and iterate from there.