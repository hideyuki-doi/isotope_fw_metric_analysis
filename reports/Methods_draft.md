# Methods (draft v1)

*Values and settings are taken from the analysis scripts (`R/01`–`R/10` in the GitHub repository). Items marked [CHECK] need confirmation before submission.*

## Data

We used the IsoFresh database [CITE], which compiles stable carbon (δ13C) and nitrogen (δ15N) isotope values for fishes and other freshwater taxa (invertebrates, primary producers and basal resources) from 1,001 food webs worldwide, together with site coordinates, ecosystem type and source publication. Taxa were identified by scientific name, or by common name where the scientific name was missing; records lacking either isotope value were removed. For each food web, isotope values of the same taxon were averaged, so that each taxon was represented by one node. Food webs with fewer than three taxa or without coordinates were excluded, leaving 981 food webs. Because network metrics are unstable in very small webs, all analyses reported here use the 606 food webs with at least 10 taxa (S = 10–96; 275 lentic and 331 lotic; 220 source studies in six regions). One web classified as both lentic and lotic was treated as lentic.

## Isotope similarity networks

Within each food web, we computed the Euclidean distance between every pair of taxa in δ13C–δ15N space and converted it to a similarity, s_ij = 1 − d_ij / max(d), where max(d) is the largest pairwise distance in that web. Pairs were linked when their similarity exceeded a single absolute threshold applied to all webs (0.662, the mean of all positive pairwise similarities pooled across webs). We used one global threshold rather than a web-specific (proportional) threshold because proportional thresholds fix link density by construction and would mask differences in connectance among webs. The resulting networks are undirected and weighted by similarity. Links therefore represent overlap in isotopic niche rather than observed feeding interactions.

## Network metrics

- **Connectance**: C = L / [S(S − 1)/2], where L is the number of links and S the number of taxa.
- **Mean link strength**: the mean, across taxa, of the summed similarity of each taxon's links (weighted degree).
- **Mean closeness**: mean normalised closeness centrality, with link lengths defined as 1 − s_ij.
- **Nestedness**: weighted NODF (vegan::nestednodf, weighted = TRUE).
- **Modularity (Q)**: estimated with simulated annealing in rnetcarto, with the number of iterations scaled to web size (iteration factor = max(10, min(100, 2000/S))).

Nestedness and modularity depend strongly on network size and link density. We therefore standardised both against a null model that randomly permuted the pairwise similarities among node pairs (preserving symmetry and the distribution of similarities) before applying the same threshold. Standardised values were computed as z = (observed − mean_null) / SD_null (z-NODF and z-Q). We used 199 null networks for webs with S ≤ 30, 100 for 30 < S ≤ 50 and 50 for S > 50 to keep computation tractable.

## Climate variables and covariates

- **Latitude**: absolute latitude of each site, used as the integrative measure of the climatic gradient.
- **Mean annual temperature**: WorldClim 2.1 BIO1 at 10 arc-minute resolution (1970–2000 normals; Fick & Hijmans 2017), extracted at each site's coordinates.
- **Temperature seasonality**: WorldClim 2.1 BIO4 divided by 100, i.e. the standard deviation of the 12 monthly mean temperatures (°C). It measures the difference between summer and winter temperatures and ranged from 0.4 to 15.9 °C (median 4.6 °C) across the webs analysed.
- **Taxon richness**: the number of isotopically characterised taxa per web (log10-transformed).
- **Prey niche area**: the convex-hull area in δ13C–δ15N space of all non-fish taxa (Layman et al. 2007).
- **Sampling effort**: IsoFresh reports, for each taxon mean, the number of samples analysed (`N_Sample`). We summed these per web and log10-transformed the total. This is usually the number of individuals but equals the number of composite samples where individuals were pooled. Values were missing for 42 webs (7%), which were excluded from analyses that included sampling effort.

Lake area (HydroLAKES; lentic webs only), the number of non-fish taxa, fish body-size dispersion and fish omnivory indices (FishBase, via rfishbase) were also compiled and analysed as exploratory covariates (Fig. S3).

Because climatic normals do not match the sampling dates of individual studies, we treat latitude as the primary descriptor of the gradient and use temperature and seasonality to interpret it.

## Statistical analyses

**Climate effects.** Each network metric was modelled as a function of one climate variable at a time using linear mixed models with nested random intercepts for region and source study, (1 | Region/Source_DOI), fitted with lme4 and lmerTest (Satterthwaite degrees of freedom). Responses and predictors were z-scored, so coefficients are standardised effects. Models were fitted to all webs and separately to lentic and lotic webs. P values were adjusted for multiple testing with the Benjamini–Hochberg false discovery rate across all tests in the model set.

**Robustness.** We repeated the climate models (i) with log10 sampling effort as an additional covariate and (ii) restricted to webs with at least the median sampling effort (n = 282). We tested sensitivity to the link threshold by recomputing connectance at thresholds between the 30th and 80th percentiles of the pooled similarity distribution. We tested whether temperature effects were modified by seasonality or by fish body-size dispersion by adding the corresponding interaction terms.

**Structural equation models.** To separate direct climate effects from effects mediated by community composition, we fitted piecewise structural equation models (piecewiseSEM; Lefcheck 2016) in which each component was a linear mixed model with the same random-effect structure:

- richness ~ climate + effort
- prey niche area ~ climate + richness + effort
- connectance, z-NODF and mean link strength ~ climate + richness + prey niche area + effort

Separate models were fitted with latitude, mean annual temperature and temperature seasonality as the climate variable, each on the same 474 webs with complete data. All variables were z-scored, so path coefficients are standardised. Indirect effects were computed as products of path coefficients along each route (climate → richness → response; climate → prey niche area → response; climate → richness → prey niche area → response), and total effects as the sum of direct and indirect effects. The same structure with Q and z-Q as responses was used to decompose modularity.

**Bootstrap mediation.** To obtain confidence intervals for indirect effects that respect the non-independence of webs from the same study, we used a study-level cluster bootstrap (999 resamples of Source_DOI with replacement). In each resample, the a path (mediator ~ climate + effort) and b path (response ~ climate + mediator + effort) were refitted as mixed models, and the indirect effect a × b was recorded. Indirect effects whose 95% percentile interval excluded zero were considered significant.

**Comparison with the temperature–connectance model.** To compare our results with Petchey, Brose & Rall (2010), we estimated the raw change in connectance per °C of mean annual temperature (mixed model with sampling effort as a covariate), both overall and with richness and prey niche area held constant. In their model, the sign of the temperature–connectance relationship depends on whether the activation energy of attack rate (E_A) or of handling time (E_H) is larger in magnitude; warming increases connectance when |E_H| > E_A. We derived the predicted sign from the activation energies reported by Rall et al. (2012) for freshwater fishes (E_A = 0.31 eV, E_H = −0.46 eV) and other groups. We also estimated the apparent activation energy of richness and connectance from mixed models of ln(Y) against inverse temperature, 1/kT (k = 8.617 × 10⁻⁵ eV K⁻¹), centred at 15 °C.

**Fish trophic position.** As an index of omnivory, we computed fish trophic height as (mean fish δ15N − mean δ15N of non-fish taxa) / 3.4 for webs with at least two fish and three non-fish taxa (n = 473 with complete data). Its relationship with climate and connectance was tested with mixed models, and it was added as a mediator of connectance in an extended SEM.

**Isotopic niche geometry.** To ask whether richer webs fill a larger isotopic space (expansion) or pack taxa more tightly (packing), we related the total convex-hull area and mean nearest-neighbour distance of all taxa to richness on log–log scales. A packing ratio was computed as the observed mean nearest-neighbour distance divided by its expectation when non-hull taxa were placed uniformly at random within the observed convex hull (hull vertices fixed; 199 randomisations).

## Software and reproducibility

Network metrics, the main models and the sampling-effort analyses were run in R 4.5.1 (lme4 2.0.6, lmerTest 3.2.1, piecewiseSEM 2.3.1, rnetcarto 0.2.6, igraph 2.2.1, vegan 2.7.2, terra 1.9.46, geodata 0.6.9, rfishbase 5.0.1). The sampling-effort analyses were also run in R 4.4.2 and gave identical results to four decimal places. The comparison with Petchey et al. (2010), fish trophic height and niche-geometry analyses were run in R 4.4.2. All code is available at https://github.com/hideyuki-doi/isotope_fw_metric_analysis [CHECK: archive a release on Zenodo for a DOI].

## References

- Fick SE, Hijmans RJ (2017) WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *Int J Climatol* 37: 4302–4315.
- Layman CA, Arrington DA, Montaña CG, Post DM (2007) Can stable isotope ratios provide for community-wide measures of trophic structure? *Ecology* 88: 42–48.
- Lefcheck JS (2016) piecewiseSEM: piecewise structural equation modelling in R for ecology, evolution, and systematics. *Methods Ecol Evol* 7: 573–579.
- Petchey OL, Brose U, Rall BC (2010) Predicting the effects of temperature on food web connectance. *Phil Trans R Soc B* 365: 2081–2091.
- Rall BC et al. (2012) Universal temperature and body-mass scaling of feeding rates. *Phil Trans R Soc B* 367: 2923–2934.
- [CITE] IsoFresh database reference.
