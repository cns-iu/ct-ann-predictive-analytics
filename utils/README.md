# Utilities

Reusable code has been modularized into functions:

1. `asctb_util_functions.py` - Data fetching and wrangling functionality using the [ASCT+B API](https://mmpyikxkcp.us-east-2.awsapprunner.com/#/).
2. `azimuth_util_functions.py` - Data fetching and wrangling functionality for the Azimuth-website's backend [repository](https://github.com/satijalab/azimuth_website/tree/master/static/csv).
3. `crosswalk_util_functions.py` - Data wrangling functionality to parse the local "crosswalk" file present in the "/asctb_mapper/" folder.
4. `utility_functions.py` - General reusable functions for summarizing the Cell-by-Gene query data, fetching top-k and bottom-k genes in the query data, analyzing pairwise predicted-annotations, etc.
5. `plotting.py` - Custom visualizations/deliverables that could be made reusable. Numerous plots in the 3 driver notebookes were not easy to functionalize.

Some functionality was already a part of the [`asctb-ct-label-mapper`](https://github.com/hubmapconsortium/asctb-ct-label-mapper/tree/main) package.