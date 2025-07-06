:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

attribute

.. auto{{ objtype }}:: {{ fullname | replace("nudca.", "nudca::") }}

{# In the fullname (e.g. `nudca.ma.MaskedArray.methodname`), the module name
is ambiguous. Using a `::` separator (e.g. `nudca::ma.MaskedArray.methodname`)
specifies `nudca` as the module name. #}