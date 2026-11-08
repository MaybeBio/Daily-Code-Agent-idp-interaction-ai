# Training Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/build_model.py](src/build_model.py)
- [src/dataset_loaders.py](src/dataset_loaders.py)
- [src/hparams_search.py](src/hparams_search.py)
- [src/loss.py](src/loss.py)
- [src/metrics.py](src/metrics.py)
- [src/utils.py](src/utils.py)

</details>



## Purpose and Scope

This document describes the training pipeline for Disobind models, focusing on the `Trainer` class and its orchestration of the model training process. The training pipeline handles optimization, learning rate scheduling, loss calculation, metrics tracking, and model calibration. For details about the model architecture being trained, see [Epsilon_3 Network Design](4.1). For information about model configuration parameters, see [Model Configuration Files](4.2). For practical guidance on training your own models, see [Training Your Own Models](4.5).

---

## Trainer Class Architecture

The `Trainer` class ([src/build_model.py:25-619]()) is the central orchestrator for model training. It inherits from `nn.Module` and manages the complete training workflow including forward passes, backpropagation, optimization, validation, and calibration.

### Core Components

```mermaid
graph TB
    subgraph "Trainer Class Initialization"
        Config["Training Configuration<br/>from YAML"]
        TrainerInit["Trainer.__init__()<br/>src/build_model.py:26-48"]
        
        Config --> TrainerInit
        
        TrainerInit --> Objective["objective: task type<br/>interaction/interface"]
        TrainerInit --> Optimizer["optim: AdamW/Adam/SGD"]
        TrainerInit --> Scheduler["scheduler_config<br/>exp/linear/multistep"]
        TrainerInit --> Loss["loss_func<br/>se_loss/bce"]
        TrainerInit --> Metrics["num_metrics: 7<br/>global averaging"]
        TrainerInit --> Calibration["method: beta/platt/iso"]
    end
    
    subgraph "Training Components"
        TrainerInit --> Model["self.model1<br/>Epsilon_3 instance"]
        TrainerInit --> Opt["self.optimizer1<br/>torch optimizer"]
        TrainerInit --> Sched["self.scheduler1<br/>LR scheduler"]
        TrainerInit --> CalModel["self.cal_model<br/>calibration model"]
    end
    
    subgraph "Training Loop Methods"
        Model --> Forward["forward()<br/>lines 512-619"]
        Forward --> TrainStep["training_step()<br/>lines 219-294"]
        Forward --> ValStep["validation_step()<br/>lines 342-400"]
        Forward --> TestStep["test_step()<br/>lines 404-491"]
    end
    
    subgraph "Helper Methods"
        TrainStep --> Predict["predict()<br/>lines 156-176"]
        TrainStep --> CalcLoss["calculate_loss_n_metrics()<br/>lines 179-215"]
        TrainStep --> GetInputs["get_inputs()<br/>lines 103-152"]
        TrainStep --> Calibrate["calibrate_model()<br/>lines 298-318"]
    end
```

**Sources:** [src/build_model.py:25-619]()

### Key Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `objective` | list | Task specification: `[task_type, bin_size, pool_type, ...]` |
| `loss_func` | function | Loss function from `get_loss_function()` [src/build_model.py:38]() |
| `optimizer1` | torch.optim | Optimizer instance (AdamW/Adam/SGD) [src/build_model.py:50-59]() |
| `scheduler1` | torch.optim.lr_scheduler | Learning rate scheduler [src/build_model.py:61-100]() |
| `cal_model` | sklearn model | Calibration model (LogisticRegression/IsotonicRegression/BetaCalibration) [src/build_model.py:298-318]() |
| `max_epochs` | int | Number of training epochs [src/build_model.py:43]() |
| `threshold` | float | Classification threshold (default 0.5) [src/build_model.py:44]() |
| `device` | str | 'cuda' or 'cpu' [src/build_model.py:45]() |

**Sources:** [src/build_model.py:26-48](), [src/build_model.py:50-100]()

---

## Training Loop

The training loop executes for `max_epochs` iterations, performing training, validation, and optionally testing at each epoch. The `HparamSearch` class [src/hparams_search.py:52]() typically manages the execution of this loop via the `Trainer.forward()` method [src/build_model.py:512]().

### Main Training Flow

```mermaid
graph TB
    Start["Trainer.forward() Entry<br/>src/build_model.py:512"] --> Initialize["Initialize<br/>optimizer, scheduler<br/>lines 541-543"]
    
    Initialize --> EpochLoop{"For each<br/>epoch"}
    
    EpochLoop -->|epoch < max_epochs| TrainMode["model.train()<br/>line 548"]
    TrainMode --> TrainStep["training_step(train_set, epoch)<br/>lines 219-294"]
    
    TrainStep --> EvalMode["model.eval()<br/>line 551"]
    EvalMode --> ValStep["validation_step(dev_set, epoch)<br/>lines 342-400"]
    
    ValStep --> RecordMetrics["Record epoch metrics<br/>lines 558-567"]
    RecordMetrics --> PrintMetrics["Print Loss, Recall,<br/>Precision, F1<br/>lines 571-576"]
    
    PrintMetrics --> EpochLoop
    
    EpochLoop -->|epoch == max_epochs| TestFinal["test_step(test_set)<br/>lines 585"]
    TestFinal --> SaveResults["Return model,<br/>cal_model, metrics<br/>line 618"]
    
    SaveResults --> End["Training Complete"]
```

**Sources:** [src/build_model.py:512-619](), [src/hparams_search.py:173-220]()

### Training Step Details

The `training_step()` method processes all mini-batches in the training set provided by the `DatasetLoader` [src/dataset_loaders.py:69]():

```mermaid
graph TB
    StartTrain["training_step(train_set, epoch)<br/>line 219"] --> InitBatch["Initialize batch_dict<br/>lines 236-237"]
    
    InitBatch --> BatchLoop{"For each<br/>batch in<br/>train_set"}
    
    BatchLoop --> ZeroGrad["optimizer.zero_grad()<br/>line 244"]
    ZeroGrad --> Predict["predict(batch, train=True)<br/>line 246"]
    
    Predict --> CalcLoss["calculate_loss_n_metrics()<br/>line 248"]
    CalcLoss --> RecordBatch["Record batch metrics<br/>lines 250-257"]
    
    RecordBatch --> Backward["loss.backward()<br/>line 260"]
    Backward --> ClipGrad{"max_norm<br/>!= None?"}
    
    ClipGrad -->|Yes| Clip["torch.nn.utils.clip_grad_norm()<br/>line 273"]
    ClipGrad -->|No| OptStep["optimizer.step()<br/>line 276"]
    Clip --> OptStep
    
    OptStep --> SWACheck{"Using SWA<br/>scheduler?"}
    SWACheck -->|Yes| UpdateSWA["Update SWA model<br/>lines 278-280"]
    SWACheck -->|No| BatchLoop
    UpdateSWA --> BatchLoop
    
    BatchLoop -->|Last epoch| CalibrateFinal["calibrate_model()<br/>lines 287-291"]
    BatchLoop -->|Not last| SchedulerStep["scheduler.step()<br/>lines 282-283"]
    
    SchedulerStep --> ReturnMetrics["Return batch_dict<br/>line 294"]
    CalibrateFinal --> ReturnMetrics
```

**Sources:** [src/build_model.py:219-294](), [src/dataset_loaders.py:106-114]()

---

## Input Preparation and Forward Pass

Before each forward pass, inputs must be prepared according to the training objective.

### Input Processing Pipeline

```mermaid
graph TB
    GetInputs["Trainer.get_inputs(batch)<br/>src/build_model.py:103"] --> SplitBatch["Split batch into<br/>target, p1_emb, p2_emb<br/>lines 125-136"]
    
    SplitBatch --> PrepInput["prepare_input()<br/>src/utils.py lines 92-205"]
    
    PrepInput --> CheckObj{"Objective<br/>type?"}
    
    CheckObj -->|interaction| MaxPool2D["Apply MaxPool2d for<br/>coarse-graining<br/>lines 144-149"]
    CheckObj -->|interface| InterfaceProj["Project to interface<br/>representation<br/>lines 169-193"]
    
    MaxPool2D --> ToDevice["Move to device<br/>lines 146-150"]
    InterfaceProj --> ToDevice
    
    ToDevice --> ReturnTensors["Return:<br/>p1_emb, p2_emb,<br/>target, masks"]
    
    ReturnTensors --> PredictForward["Trainer.predict(batch)<br/>src/build_model.py:156"]
    PredictForward --> ModelForward["model(p1_emb, p2_emb,<br/>interaction_mask)<br/>line 173"]
    
    ModelForward --> ReturnPreds["Return predictions,<br/>target, mask<br/>line 176"]
```

**Sources:** [src/build_model.py:103-176](), [src/utils.py:92-205]()

### Coarse-Graining Process

For binned (coarse-grained) predictions, the input preparation in `prepare_input` [src/utils.py:92]() applies pooling operations:

| Objective | CG Level | Pooling Type | Implementation |
|-----------|----------|--------------|----------------|
| `interaction_bin` | > 1 | `MaxPool2d` | `nn.MaxPool2d(kernel_size=bin_size, stride=bin_size)` [src/utils.py:144]() |
| `interface_bin` | > 1 | `MaxPool2d` | Applied to projected interface mask [src/utils.py:153]() |
| Embedding Averaging | > 1 | `AvgPool1d` | `nn.AvgPool1d` if `bin_input` is True [src/utils.py:136]() |

**Sources:** [src/utils.py:129-156]()

---

## Optimization Strategy

### Optimizer Configuration

The `Trainer` supports three optimizers initialized in the `optimizer()` method [src/build_model.py:50]():

```mermaid
graph LR
    OptimizerMethod["Trainer.optimizer()<br/>src/build_model.py:50"] --> CheckType{"optimizer<br/>type?"}
    
    CheckType -->|Adam| AdamOpt["torch.optim.Adam<br/>lr, weight_decay, amsgrad<br/>line 52"]
    CheckType -->|SGD| SGDOpt["torch.optim.SGD<br/>lr, weight_decay, momentum=0.9<br/>line 55"]
    CheckType -->|AdamW| AdamWOpt["torch.optim.AdamW<br/>lr, weight_decay, amsgrad<br/>line 58"]
    
    AdamOpt --> Return["Return optimizer"]
    SGDOpt --> Return
    AdamWOpt --> Return
```

**Sources:** [src/build_model.py:50-59]()

### Learning Rate Schedulers

The `scheduler()` method ([src/build_model.py:61-100]()) configures various learning rate schedules based on the `scheduler_config` [src/build_model.py:33]():

| Scheduler | Type | Key Parameters |
|-----------|------|----------------|
| `swa` | Stochastic Weight Averaging | `swa_lr`, `swa_start` [src/build_model.py:64-68]() |
| `cycliclr` | Cyclic Learning Rate | `base_lr`, `max_lr`, `step_size_up` [src/build_model.py:69-75]() |
| `multistep` | MultiStep Decay | `milestones`, `gamma` [src/build_model.py:76-80]() |
| `exp` | Exponential Decay | `gamma` [src/build_model.py:82-85]() |
| `linear` | Linear Decay | `start_factor`, `end_factor`, `total_iters` [src/build_model.py:87-94]() |

**Sources:** [src/build_model.py:61-100]()

---

## Loss Functions

Disobind implements several loss functions in `src/loss.py` to handle the specific requirements of protein interaction prediction.

| Loss Class | Implementation | Use Case |
|------------|----------------|----------|
| `BCE` | `nn.BCELoss` | Standard binary cross entropy [src/loss.py:55]() |
| `BCEwithLogits` | `nn.BCEWithLogitsLoss` | Weighted BCE with logit input [src/loss.py:104]() |
| `FocalLoss` | Custom implementation | Handles class imbalance by focusing on hard examples [src/loss.py:135]() |
| `SingularityEnhancedLoss` | Custom implementation | Specifically designed for contact prediction [src/loss.py:162]() |
| `InterfaceLoss` | Wrapper | Combines losses for dual interface predictions [src/loss.py:189]() |

The `get_loss_function()` helper in `src/loss.py` (referenced in [src/build_model.py:38]()) instantiates these classes based on the configuration.

**Sources:** [src/loss.py:1-205](), [src/build_model.py:38]()

---

## Validation and Testing

### Validation Step

The `validation_step()` method [src/build_model.py:342]() evaluates the model on the development set without gradient updates. It also collects predictions for calibration during the final epoch [src/build_model.py:371-379]().

### Test Step

The `test_step()` method [src/build_model.py:404]() generates final predictions on the test set and produces reliability diagrams using `plot_reliabity_diagram` [src/utils.py:21]().

```mermaid
graph TB
    TestStart["test_step(test_set, file_name)<br/>line 404"] --> NoGradTest["with torch.no_grad():<br/>line 428"]
    
    NoGradTest --> TestLoop{"For each<br/>batch"}
    
    TestLoop --> PredTest["predict(batch)<br/>line 432"]
    PredTest --> GetCalPreds["get_calibrated_preds()<br/>lines 447-454"]
    
    GetCalPreds --> CalcMetricsTest["calculate_loss_n_metrics()<br/>line 460"]
    CalcMetricsTest --> RecordMetricsTest["Record metrics<br/>lines 463-469"]
    
    RecordMetricsTest --> TestLoop
    
    TestLoop -->|Complete| AvgMetrics["Average metrics<br/>line 474"]
    AvgMetrics --> CheckMethod{"Calibration<br/>method?"}
    
    CheckMethod -->|None| SkipPlot["Skip reliability plot<br/>line 481"]
    CheckMethod -->|temp/beta/platt/iso| PlotReliability["plot_reliabity_diagram()<br/>src/utils.py:21"]
    
    SkipPlot --> PrintResults["Print Recall, Precision, F1<br/>lines 488-489"]
    PlotReliability --> PrintResults
    
    PrintResults --> Return["Return test_met,<br/>test_pred, test_target<br/>line 491"]
```

**Sources:** [src/build_model.py:404-491](), [src/utils.py:21-62]()

---

## Model Calibration

Calibration improves the reliability of predicted probabilities. The `Trainer` supports multiple calibration methods via `calibrate_model` [src/build_model.py:298]().

### Calibration Methods Comparison

| Method | Implementation | Library |
|--------|----------------|---------|
| `platt` | `LogisticRegression` | `sklearn.linear_model` [src/build_model.py:302]() |
| `iso` | `IsotonicRegression` | `sklearn.isotonic` [src/build_model.py:307]() |
| `beta-abm` | `BetaCalibration` | `betacal` [src/build_model.py:311]() |
| `temp` | Temperature Scaling | Internal [src/build_model.py:318]() |

Calibrated predictions are retrieved using `get_calibrated_preds()` [src/build_model.py:322](), which applies the fitted `self.cal_model`.

**Sources:** [src/build_model.py:298-338]()

---