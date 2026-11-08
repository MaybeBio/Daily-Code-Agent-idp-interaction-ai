# GUI Architecture and Components

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [GMXMMPBSA/analyzer/chartsettings.py](GMXMMPBSA/analyzer/chartsettings.py)
- [GMXMMPBSA/analyzer/customitem.py](GMXMMPBSA/analyzer/customitem.py)
- [GMXMMPBSA/analyzer/dialogs.py](GMXMMPBSA/analyzer/dialogs.py)
- [GMXMMPBSA/analyzer/gui.py](GMXMMPBSA/analyzer/gui.py)
- [GMXMMPBSA/analyzer/items_delegate.py](GMXMMPBSA/analyzer/items_delegate.py)
- [GMXMMPBSA/analyzer/parametertree/Parameter.py](GMXMMPBSA/analyzer/parametertree/Parameter.py)
- [GMXMMPBSA/analyzer/parametertree/ParameterItem.py](GMXMMPBSA/analyzer/parametertree/ParameterItem.py)
- [GMXMMPBSA/analyzer/parametertree/ParameterTree.py](GMXMMPBSA/analyzer/parametertree/ParameterTree.py)
- [GMXMMPBSA/analyzer/parametertree/__init__.py](GMXMMPBSA/analyzer/parametertree/__init__.py)
- [GMXMMPBSA/analyzer/parametertree/parameterTypes.py](GMXMMPBSA/analyzer/parametertree/parameterTypes.py)
- [GMXMMPBSA/analyzer/plots.py](GMXMMPBSA/analyzer/plots.py)
- [GMXMMPBSA/analyzer/utils.py](GMXMMPBSA/analyzer/utils.py)

</details>



## Purpose and Scope

This page describes the internal architecture of the **gmx_MMPBSA_ana** graphical analyzer application, focusing on its GUI components, window management, and configuration system. The gmx_MMPBSA_ana application provides an interactive GUI for analyzing MM/PB(GB)SA calculation results through charts, tables, and 3D molecular visualization. It implements a Multiple Document Interface (MDI) architecture with dockable configuration panels.

---

## Main Window Architecture

### GMX_MMPBSA_ANA Class

The core application is implemented in the `GMX_MMPBSA_ANA` class, which inherits from `QMainWindow` and orchestrates all GUI components.

**Main Window Component Hierarchy**

```mermaid
graph TB
    MainWindow["GMX_MMPBSA_ANA<br/>(QMainWindow)"]
    
    MDI["QMdiArea<br/>Central Widget<br/>mdi"]
    
    TreeDock["QDockWidget<br/>treeDockWidget<br/>Right Area"]
    CorrDock["QDockWidget<br/>correlation_DockWidget<br/>Right Area"]
    OptDock["QDockWidget<br/>optionDockWidget<br/>Left Area"]
    
    TreeWidget["QTreeWidget<br/>treeWidget<br/>System hierarchy"]
    CorrContainer["QWidget<br/>corr_container_widget"]
    OptWidget["QWidget<br/>optionWidget"]
    
    MenuBar["QMenuBar<br/>menubar"]
    StatusBar["QStatusBar<br/>statusbar"]
    
    MainWindow --> MDI
    MainWindow --> TreeDock
    MainWindow --> CorrDock
    MainWindow --> OptDock
    MainWindow --> MenuBar
    MainWindow --> StatusBar
    
    TreeDock --> TreeWidget
    CorrDock --> CorrContainer
    OptDock --> OptWidget
    
    subgraph "MDI Area"
        SubWin1["Chart SubWindows<br/>LineChart, BarChart, etc."]
        SubWin2["Table SubWindows<br/>Tables, OutputFiles"]
    end
    
    MDI --> SubWin1
    MDI --> SubWin2
```

**Key Attributes:**

| Attribute | Type | Purpose |
|-----------|------|---------|
| `mdi` | `QMdiArea` | Central widget hosting multiple chart/table subwindows [GMXMMPBSA/analyzer/gui.py:61]() |
| `treeDockWidget` | `QDockWidget` | Right dock containing system/data tree [GMXMMPBSA/analyzer/gui.py:65]() |
| `correlation_DockWidget` | `QDockWidget` | Right dock for correlation analysis [GMXMMPBSA/analyzer/gui.py:68]() |
| `optionDockWidget` | `QDockWidget` | Left dock for chart options [GMXMMPBSA/analyzer/gui.py:73]() |
| `systems` | `dict` | Stores loaded system data and configuration [GMXMMPBSA/analyzer/gui.py:53]() |
| `items_counter` | `dict` | Tracks chart instances [GMXMMPBSA/analyzer/gui.py:59]() |

Sources: [GMXMMPBSA/analyzer/gui.py:44-99]()

---

## Dockable Panel System

### Data Tree Dock

The data tree dock displays the hierarchical structure of loaded systems and their energy components. It uses a `QTreeWidget` with custom items (`CustomItem`) that provide chart generation buttons.

**Data Tree Structure**

```mermaid
graph TB
    TreeWidget["QTreeWidget<br/>treeWidget"]
    
    SysItem1["System 1<br/>CustomItem"]
    SysItem2["System 2<br/>CustomItem"]
    
    EnergyItem["Energy<br/>CustomItem"]
    DecompItem["Decomposition<br/>CustomItem"]
    EntropyItem["Entropy<br/>CustomItem"]
    
    GBItem["GB<br/>CustomItem"]
    PBItem["PB<br/>CustomItem"]
    
    Toolbar["QToolBar<br/>Chart action buttons"]
    LineBtn["Line Chart Button"]
    BarBtn["Bar Chart Button"]
    HeatmapBtn["Heatmap Button"]
    VisBtn["PyMOL Visualization"]
    
    TreeWidget --> SysItem1
    TreeWidget --> SysItem2
    
    SysItem1 --> EnergyItem
    SysItem1 --> DecompItem
    SysItem1 --> EntropyItem
    
    EnergyItem --> GBItem
    EnergyItem --> PBItem
    
    GBItem --> Toolbar
    Toolbar --> LineBtn
    Toolbar --> BarBtn
    Toolbar --> HeatmapBtn
    Toolbar --> VisBtn
```

Each `CustomItem` in the tree has:
- `buttons`: Tuple defining which chart types are available (1=Line, 2=Bar, 3=Heatmap, 4=Visualization) [GMXMMPBSA/analyzer/customitem.py:114]()
- `keys_path`: Tuple path to locate data in `systems[index]['items_data']` [GMXMMPBSA/analyzer/customitem.py:118]()
- `system_index`: Reference to parent system [GMXMMPBSA/analyzer/customitem.py:109]()
- Chart subwindow references: `lp_subw`, `bp_subw`, `hmp_subw` [GMXMMPBSA/analyzer/customitem.py:120-122]()

Sources: [GMXMMPBSA/analyzer/customitem.py:104-164](), [GMXMMPBSA/analyzer/gui.py:76-99]()

### Correlation Dock

The correlation dock enables multi-system correlation analysis, displaying correlation coefficients and regression plots.

**Correlation Dock Components**

```mermaid
graph TB
    CorrDock["correlation_DockWidget"]
    
    Container["corr_container_widget<br/>QWidget"]
    
    CorrTable["correlation_tableWidget<br/>QTableWidget<br/>Model selection"]
    DataTable["data_table_widget<br/>QTableWidget<br/>System data display"]
    
    TabWidget["corr_options_widget<br/>QTabWidget"]
    
    ChartTab["Charts Tab<br/>corr_parm_tree_w<br/>ParameterTree"]
    SelectTab["Selections Tab<br/>corr_sys_sel_w<br/>System selector"]
    
    UpdateBtn["update_corr_btn<br/>QPushButton"]
    
    CorrDock --> Container
    Container --> CorrTable
    Container --> DataTable
    Container --> TabWidget
    
    TabWidget --> ChartTab
    TabWidget --> SelectTab
    TabWidget --> UpdateBtn
```

The correlation system uses `CustomCorrItem` objects that generate `RegChart` visualizations when selected [GMXMMPBSA/analyzer/customitem.py:71]().

Sources: [GMXMMPBSA/analyzer/gui.py:100-168](), [GMXMMPBSA/analyzer/customitem.py:52-102]()

### Options Dock

The options dock provides real-time chart customization through a tabbed interface with ParameterTree widgets.

**Options Dock Structure**

```mermaid
graph TB
    OptDock["optionDockWidget"]
    
    OptWidget["optionWidget<br/>QWidget"]
    
    SelGroup["Selection RadioButtons<br/>Current/All Systems"]
    
    TabWidget["optionWidget_c<br/>QTabWidget"]
    
    ChartsTab["Charts Tab"]
    FramesTab["Frames Tab"]
    
    ParmTree["parm_tree_w<br/>ParameterTree<br/>Dynamic configuration"]
    
    EFrames["Energy Frame Range<br/>Start/End/Interval"]
    NMFrames["NMODE Frame Range"]
    IESegment["IE Segment"]
    
    UpdateBtn["update_btn<br/>Apply changes"]
    
    OptDock --> OptWidget
    OptWidget --> SelGroup
    OptWidget --> TabWidget
    
    TabWidget --> ChartsTab
    TabWidget --> FramesTab
    TabWidget --> UpdateBtn
    
    ChartsTab --> ParmTree
    FramesTab --> EFrames
    FramesTab --> NMFrames
    FramesTab --> IESegment
```

Sources: [GMXMMPBSA/analyzer/gui.py:400-592]()

---

## Chart System Architecture

### Chart Class Hierarchy

All chart types inherit from `ChartsBase`, which provides common functionality for MDI subwindows and matplotlib integration [GMXMMPBSA/analyzer/plots.py:122]().

**Chart Type Hierarchy**

```mermaid
graph TB
    Base["ChartsBase<br/>(QMdiSubWindow)"]
    
    Line["LineChart<br/>Time series energy plots"]
    Bar["BarChart<br/>Energy term comparisons"]
    Heatmap["HeatmapChart<br/>Per-residue/Per-wise"]
    Reg["RegChart<br/>Correlation analysis"]
    Tables["Tables<br/>Tabular data display"]
    Output["OutputFiles<br/>Text file viewer"]
    
    Base --> Line
    Base --> Bar
    Base --> Heatmap
    Base --> Reg
    Base --> Tables
    Base --> Output
    
    MHeatmap["MHeatmap<br/>Custom heatmap renderer"]
    Heatmap --> MHeatmap
```

**ChartsBase Common Interface**

```mermaid
graph LR
    ChartsBase["ChartsBase"]
    
    FigCanvas["FigureCanvas<br/>figure_canvas<br/>Matplotlib backend"]
    NavBar["NavigationToolbar<br/>mpl_toolbar<br/>Custom save/zoom"]
    MainWidget["QMainWindow<br/>mainwidgetmdi<br/>Toolbar container"]
    
    Button["QToolButton<br/>button<br/>Toggle reference"]
    ItemParent["CustomItem<br/>item_parent<br/>Tree reference"]
    
    ChartsBase --> FigCanvas
    ChartsBase --> NavBar
    ChartsBase --> MainWidget
    ChartsBase --> Button
    ChartsBase --> ItemParent
    
    FigCanvas --> MainWidget
    NavBar --> MainWidget
```

Key methods:
- `set_cw(fig)`: Initializes figure canvas and toolbar [GMXMMPBSA/analyzer/plots.py:137-150]()
- `setup_text(ax, options, ...)`: Configures axis labels and fonts [GMXMMPBSA/analyzer/plots.py:164-187]()
- `draw()`: Renders the figure [GMXMMPBSA/analyzer/plots.py:159-162]()

Sources: [GMXMMPBSA/analyzer/plots.py:122-191]()

### Chart Creation Flow

**Data-to-Chart Pipeline**

```mermaid
graph TB
    User["User clicks chart button"]
    CustomItem["CustomItem.plotting_*()<br/>line/bar/heatmap"]
    
    CheckState["Check for existing subwindow"]
    CheckChanges["Check chart_options.changes"]
    
    LoadData["Load data from<br/>systems[idx]['items_data'][keys_path]"]
    
    CreateChart["Instantiate Chart class<br/>LineChart/BarChart/HeatmapChart"]
    
    AddMDI["mdi.addSubWindow(chart)"]
    ShowChart["chart.show()"]
    
    UpdateChart["chart.update_config(options)"]
    
    User --> CustomItem
    CustomItem --> CheckState
    CheckState -->|"Not exists or<br/>action==3"| LoadData
    CheckState -->|"Exists and<br/>action==1"| UpdateChart
    
    LoadData --> CreateChart
    CreateChart --> AddMDI
    AddMDI --> ShowChart
```

The `changes` dict tracks what needs updating: `action_type=1` for update existing, `action_type=3` for re-plot [GMXMMPBSA/analyzer/customitem.py:381-454]().

---

## ParameterTree Configuration System

The ParameterTree provides type-safe, hierarchical configuration with immediate visual feedback.

### ParameterTree Architecture

```mermaid
graph TB
    Tree["ParameterTree<br/>(QTreeWidget)"]
    
    Param["Parameter<br/>(QObject)<br/>Base parameter class"]
    
    Item["ParameterItem<br/>(QTreeWidgetItem)<br/>Display delegate"]
    
    WidgetItem["WidgetParameterItem<br/>Editor widgets"]
    
    IntType["SimpleParameter<br/>type='int'"]
    FloatType["SimpleParameter<br/>type='float'"]
    BoolType["SimpleParameter<br/>type='bool'"]
    ListType["ListParameter<br/>type='list'"]
    ColorType["ColorParameter<br/>type='color'"]
    
    Tree --> Param
    Param --> Item
    Item --> WidgetItem
    
    WidgetItem --> IntType
    WidgetItem --> FloatType
    WidgetItem --> BoolType
    WidgetItem --> ListType
    WidgetItem --> ColorType
```

### Parameter Registration and Types

The system uses a global registry to map type names to Parameter classes:

```mermaid
graph LR
    Registry["PARAM_TYPES dict"]
    
    Int["'int' → SimpleParameter"]
    Float["'float' → SimpleParameter"]
    Bool["'bool' → SimpleParameter"]
    Str["'str' → SimpleParameter"]
    List["'list' → ListParameter"]
    Color["'color' → ColorParameter"]
    Group["'group' → GroupParameter"]
    Action["'action' → ActionParameter"]
    
    Registry --> Int
    Registry --> Float
    Registry --> Bool
    Registry --> Str
    Registry --> List
    Registry --> Color
    Registry --> Group
    Registry --> Action
```

Parameter creation: `Parameter.create(**opts)` looks up `opts['type']` in `PARAM_TYPES` and instantiates the appropriate class [GMXMMPBSA/analyzer/parametertree/Parameter.py:119-129]().

Sources: [GMXMMPBSA/analyzer/parametertree/Parameter.py:30-50](), [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:120-160]()

### ChartSettings Configuration

The `ChartSettings` class wraps a nested dictionary structure that maps to ParameterTree nodes [GMXMMPBSA/analyzer/chartsettings.py:299](). Access pattern: `options[('Line Plot', 'line-width')]` returns the value [GMXMMPBSA/analyzer/chartsettings.py:1-21]().

Sources: [GMXMMPBSA/analyzer/chartsettings.py:299-800]()

---

## Widget Customization System

### Custom Widgets for Parameters

The ParameterTree system uses specialized widgets for different parameter types:

| Parameter Type | Widget Class | Purpose |
|----------------|--------------|---------|
| `int`, `float` | `QSpinBox`, `QDoubleSpinBox` | Numeric input with bounds [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:130-132]() |
| `bool` | `QCheckBox` | Toggle options [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:133]() |
| `str` | `QLineEdit` | Text input [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:134]() |
| `color` | `ColorButton` | Color picker dialog [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:135]() |

The `ColorButton` widget provides custom color selection via `QColorDialog` [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:53]().

Sources: [GMXMMPBSA/analyzer/parametertree/parameterTypes.py:35-118]()

### NavigationToolbar Customization

The analyzer overrides matplotlib's `NavigationToolbar2QT` to control save behavior and DPI settings [GMXMMPBSA/analyzer/plots.py:60]().

```mermaid
graph TB
    NavBar["NavigationToolbar<br/>(NavigationToolbar2QT)"]
    
    Attrs["Attributes:<br/>save_format<br/>save_dpi<br/>filename"]
    
    UpdateOpts["update_options(options)<br/>Set from ChartSettings"]
    
    SaveFig["save_figure()<br/>Custom save dialog<br/>with configured format"]
    
    NavBar --> Attrs
    NavBar --> UpdateOpts
    NavBar --> SaveFig
```

Sources: [GMXMMPBSA/analyzer/plots.py:60-120]()

---

## Dialog System

### InitDialog Startup Configuration

The `InitDialog` class provides initial system configuration, including energy options and decomposition settings [GMXMMPBSA/analyzer/dialogs.py:41](). It also includes a performance slider to manage resource allocation [GMXMMPBSA/analyzer/dialogs.py:145]().

Sources: [GMXMMPBSA/analyzer/dialogs.py:36-266]()

### ProcessingProgressBar

The `ProcessingProgressBar` dialog displays progress during data recalculation using multiprocessing pools [GMXMMPBSA/analyzer/dialogs.py:578]().

Sources: [GMXMMPBSA/analyzer/dialogs.py:578-700]()

---

## Summary Table: Key Classes

| Class | Module | Purpose |
|-------|--------|---------|
| `GMX_MMPBSA_ANA` | `gui.py` | Main window, MDI orchestration [GMXMMPBSA/analyzer/gui.py:44]() |
| `ChartsBase` | `plots.py` | Base class for all chart types [GMXMMPBSA/analyzer/plots.py:122]() |
| `CustomItem` | `customitem.py` | Tree widget items with chart buttons [GMXMMPBSA/analyzer/customitem.py:104]() |
| `ParameterTree` | `parametertree/ParameterTree.py` | Configuration tree widget [GMXMMPBSA/analyzer/parametertree/ParameterTree.py:30]() |
| `Parameter` | `parametertree/Parameter.py` | Configuration node base class [GMXMMPBSA/analyzer/parametertree/Parameter.py:52]() |
| `InitDialog` | `dialogs.py` | Startup configuration dialog [GMXMMPBSA/analyzer/dialogs.py:41]() |

Sources: Multiple files cited above.

---