# 0D chemicals at sea weathering model
Simple python module able to modelize some chemicals weathering at sea surface, without any spatial input.

## Table of contents
* General info
* Technologies
* Contacts
* References

## General info
The module can be used to model HNS substance or oil weathering. It needs some input such as the temperature or the windspeed. It allows for pure chemicals and for pseudo-component chemicals by providing an object Mix and an object Component (in oil_utils.py), and some equation to compute the evolution of it (in weathering_utils.py). There is a file "exemple.py" wich allows to understand how the code is working. There is also a list of parametrization function in the files describing the weathering processes (evaporation.py, dissolution.py, photooxidation.py, emulsification.py, biodegradation.py) which can be used. All the sources for these parametrization are in the references section.

## Technologies
The project need python and numpy to work.

## Contacts

## License
Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## References
Arey, J.S., Nelson, R.K., Plata, D.L., Reddy, C.M., 2007. Disentangling Oil Weathering Using GC×GC. 2. Mass Transfer Calculations. Environ. Sci. Technol. 41, 5747–5755. https://doi.org/10.1021/es070006p

Berry, A., Dabrowski, T., Lyons, K., 2012. The oil spill model OILTRANS and its application to the Celtic Sea. Marine Pollution Bulletin 64, 2489–2501. https://doi.org/10.1016/j.marpolbul.2012.07.036

Brighton P.W.M., 1985. Evaporation from a Plane Liquid Surface into a Turbulent Boundary-Layer. Journal of Fluid Mechanics 159, 323–345. https://doi.org/10.1017/S0022112085003238

CHEMMAP technical User’s manual 6.10, 2014.

Cohen, Y., Mackay, D., Shiu, W.Y., 1980. Mass transfer rates between oil slicks and water. Can. J. Chem. Eng. 58, 569–575. https://doi.org/10.1002/cjce.5450580504

D-WAQ PART User Manual, n.d. 150.
Eley, D.D., Hey, M.J., Symonds, J.D., 1988. Emulsions of water in asphaltene-containing oils 1. Droplet size distribution and emulsification rates. Colloids and Surfaces 32, 87–101. https://doi.org/10.1016/0166-6622(88)80006-4

Fingas, M. (Ed.), 2015. Handbook of oil spill science and technology, 1. ed. ed. Wiley, Hoboken, NJ.

Fingas, M., 1995. Water-in-oil emulsion formation: A review of physics and mathematical modelling. Spill Science & Technology Bulletin 2, 55–59. https://doi.org/10.1016/1353-2561(95)94483-Z

Jones, R., Lehr, W., Simecek-Beatty, D., Reynolds, R.M., n.d. ALOHA® (Areal Locations Of Hazardous Atmospheres) 5.4.4: Technical Documentation 96.

Keramea, P., Spanoudaki, K., Zodiatis, G., Gikas, G., Sylaios, G., 2021. Oil Spill Modeling: A Critical Review on Current Trends, Perspectives, and Challenges. JMSE 9, 181. https://doi.org/10.3390/jmse9020181

Kotzakoulakis, K., George, S.C., 2018. Predicting the weathering of fuel and oil spills: A diffusion-limited evaporation model. Chemosphere 190, 442–453. https://doi.org/10.1016/j.chemosphere.2017.09.142

Lehr, W., Jones, R., Evans, M., Simecek-Beatty, D., Overstreet, R., 2002. Revisions of the ADIOS oil spill model. Environmental Modelling & Software 17, 189–197. https://doi.org/10.1016/S1364-8152(01)00064-0

Mackay, D., Matsugu, R.S., 1973. Evaporation rates of liquid hydrocarbon spills on land and water. Can. J. Chem. Eng. 51, 434–439. https://doi.org/10.1002/cjce.5450510407

Mishra, A.K., Kumar, G.S., 2015. Weathering of Oil Spill: Modeling and Analysis. Aquatic Procedia 4, 435–442. https://doi.org/10.1016/j.aqpro.2015.02.058

Nordam, T., 2020. Modelling biodegradation of crude oil components at low temperatures 4.

Shen, H.T., Yapa, P.D., Wang, D.S., Yang, X.Q., 1993. A Mathematical Model for Oil Slick Transport and Mixing in Rivers 79.

Stiver, Warren., Mackay, Donald., 1984. Evaporation rate of spills of hydrocarbons and petroleum mixtures. Environ. Sci. Technol. 18, 834–840. https://doi.org/10.1021/es00129a006
Sutton, O.G., Simpson, G.C., 1934. Wind structure and evaporation in a turbulent atmosphere. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character 146, 701–722. https://doi.org/10.1098/rspa.1934.0183
Vaz, A.C., Faillettaz, R., Paris, C.B., 2021. A Coupled Lagrangian-Earth System Model for Predicting Oil Photooxidation. Front. Mar. Sci. 8, 576747. https://doi.org/10.3389/fmars.2021.576747

Vilcáez, J., Li, L., Hubbard, S.S., 2013. A new model for the biodegradation kinetics of oil droplets: application to the Deepwater Horizon oil spill in the Gulf of Mexico. Geochem Trans 14, 4. https://doi.org/10.1186/1467-4866-14-4

Ward, C.P., Armstrong, C.J., Conmy, R.N., French-McCay, D.P., Reddy, C.M., 2018. Photochemical Oxidation of Oil Reduced the Effectiveness of Aerial Dispersants Applied in Response to the Deepwater Horizon Spill. Environ. Sci. Technol. Lett. 5, 226–231. https://doi.org/10.1021/acs.estlett.8b00084
