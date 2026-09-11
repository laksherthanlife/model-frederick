"""The mechanistic layer: ODE state between the environment and the existing product code.

Nothing is imported here on purpose. The modules of this package are being written
independently and an eager re-export would make each one's import errors the whole
package's, which is the failure mode `ystwin/__init__.py` already avoids by carrying a
version string and nothing else.

`params.py` is the only file every other module in here is expected to depend on. It
carries :class:`~ystwin.mech.params.Param` and the tag taxonomy that makes criterion (b) --
every state constrained by something measurable, or a DECLARED sweep axis -- checkable by a
test rather than by reading a docstring, and the free-scalar gate of criterion (e).
"""
