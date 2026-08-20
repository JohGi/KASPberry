# Vendored blastn2dotplots

This directory contains a modified copy of blastn2dotplots:
https://github.com/mokuno3430/blastn2dotplots

Original software is distributed under the MIT License.
See LICENSE.

KASPberry modifications:
- generation of an additional cropped, decoration-free dotplot PDF
  for embedding in the region viewer.

The upstream `paf2blastn-fmt6.pl` utility is also retained because
KASPberry uses it to convert minimap2 PAF output to the input format
expected by blastn2dotplots.
