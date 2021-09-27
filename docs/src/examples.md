# Examples

```@setup highlight
using CodeTracking, Markdown, OpticSim.Examples
mdparse(s) = Markdown.parse("```julia\n$s\n```")
```

```@setup example
using OpticSim.Examples
```

## Pluto Notebooks

The [`OpticSim`](index.html) package comes with several [`Pluto`](https://github.com/fonsp/Pluto.jl) notebooks (code snippets are coming soon) that allow the user to change and run sample code and view the results in real-time. We highly recommend for you to try these out.
The notebooks are located in the **samples** folder, and you can try them by running:

```julia
import OpticSim.NotebooksUtils as NB    # the **as** option was added in Julia v1.6

NB.run_sample("EmittersIntro.jl")
```

The **run_sample** method will **copy** the notebook to your current folder (if it does not exist) and launch Pluto to run the notebook in the browser.

## Cooke Triplet
```@example highlight
mdparse(@code_string draw_cooketriplet()) # hide
```

```@example example
sys = draw_cooketriplet("assets/cooke.png"); @show sys; nothing # hide
```

![Cooke triplet visualization](assets/cooke.png)

## Zoom Lens
```@example highlight
mdparse(@code_string draw_zoomlenses()) # hide
```

```@example example
syss = draw_zoomlenses(["assets/zoom$i.png" for i in 1:3]); @show syss[1]; nothing # hide
```

![Zoom position 1 visualization](assets/zoom1.png)
![Zoom position 2 visualization](assets/zoom2.png)
![Zoom position 3 visualization](assets/zoom3.png)

## Schmidt Cassegrain Telescope
```@example highlight
mdparse(@code_string draw_schmidtcassegraintelescope()) # hide
```

```@example example
draw_schmidtcassegraintelescope("assets/tele.png") # hide
```

![Schmidt Cassegrain Telescope visualization](assets/tele.png)

## Lens Construction
```@example highlight
mdparse(@code_string draw_lensconstruction()) # hide
```

```@example example
draw_lensconstruction("assets/qtype.png") # hide
```

![lens construction example](assets/qtype.png)

## HOEs

### Focusing
```@example highlight
mdparse(@code_string draw_HOEfocus()) # hide
```

```@example example
draw_HOEfocus("assets/hoe_f.png") # hide
```

![Focusing HOE example](assets/hoe_f.png)

### Collimating
```@example highlight
mdparse(@code_string draw_HOEcollimate()) # hide
```

```@example example
draw_HOEcollimate("assets/hoe_c.png") # hide
```

![Collimating HOE example](assets/hoe_c.png)

### Multi
```@example highlight
mdparse(@code_string draw_multiHOE()) # hide
```

```@example example
draw_multiHOE("assets/hoe_m.png") # hide
```

![Multi-HOE example](assets/hoe_m.png)

## Deterministic Raytracing
```@example highlight
mdparse(@code_string draw_stackedbeamsplitters()) # hide
```

```@example example
draw_stackedbeamsplitters(["assets/deterministic_trace_$i.png" for i in 1:3]) # hide
```

![Nondeterministic Raytrace](assets/deterministic_trace_1.png)
![Transmission only](assets/deterministic_trace_2.png)
![Reflection only](assets/deterministic_trace_3.png)

## Gestel
```@example highlight
mdparse(@code_string draw_gestel()) # hide
```

```@example
using OpticSim, OpticSim.Examples
OpticSim.NotebooksUtils.SetDocsBackend("Web")
draw_gestel()
```
