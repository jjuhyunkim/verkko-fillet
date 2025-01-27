## Preprocessing : `pp`

```{eval-rst}
.. module:: verkkofillet.pp
```

```{eval-rst}
.. currentmodule:: verkkofillet
```

Extracting gap from path, preprocessing recipes. 

### Reading and writing
```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.read_Verkko
   pp.readChr
   pp.writeFixedGaf
   pp.save_Verkko
   pp.load_Verkko
   pp.mkCNSdir
   pp.updateCNSdir_missingEdges

``` 
### Preprocessing
```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.find_intra_telo
   pp.find_reads_intra_telo
   pp.findGaps
   pp.searchNodes
   pp.highlight_nodes
   pp.fillGaps
   pp.estLoops
   pp.checkGapFilling
   pp.naming_contigs
   pp.find_multi_used_node

``` 