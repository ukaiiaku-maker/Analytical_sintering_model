# Defined closed-pore law implementation

`mechanism_mode="defined_closed_laws_port"` assembles the already-defined open renewal, Onsager stress, conservative PR/closure, binwise finite accommodation, closed renewal/GB/gas alternatives, non-densifying surface accommodation, and migration-only Zener/growth laws. It introduces no new mechanism.

Closed bins track volume, count proxy, evolving radius, shrinkability, geometry, gas pressure, accommodation maximum/used/recovered/available, and named closed removal. The default renewal radius exponent is four; three is a sensitivity. The fitted ZrO2 barrier and the 380 kJ/mol GB/surface diffusivities are unchanged. Local functions receive state, temperature, material, and geometry only.

The discussion DOCX was unavailable locally; its controlling rules were taken from the supplied correction brief. Available renewal/Onsager Python and MATLAB scripts provide equation-level traceability. No validation is claimed.
