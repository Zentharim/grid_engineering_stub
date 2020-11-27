# KEEP core/proj/reproj


class OverlapError(Exception):
	def __init__(self):
		pass
	def __str__(self):
		return "Non overlap data"
# KEEP core/proj/reproj
class ReprojError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)
