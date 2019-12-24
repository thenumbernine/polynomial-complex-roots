local class = require 'ext.class'
local math = require 'ext.math'

local cplx = class()

function cplx:init(re,im)
	if type(re) == 'table' and cplx.is(re) then
		self[1] = re[1] or 0
		self[2] = re[2] or 0
	else
		self[1] = re or 0
		self[2] = im or 0
	end
end

function cplx.from(a)
	if type(a) == 'table' then
		return cplx(a)
	elseif type(a) == 'number' then
		return cplx(a,0)
	end
	return error("couldn't convert to cplx")
end

function cplx.__unm(a)
	a = cplx.from(a)
	return cplx(-a[1], -a[2])
end

function cplx.__add(a,b)
	a = cplx.from(a)
	b = cplx.from(b)
	return cplx(a[1] + b[1], a[2] + b[2])
end

function cplx.__sub(a,b)
	a = cplx.from(a)
	b = cplx.from(b)
	return cplx(a[1] - b[1], a[2] - b[2])
end

function cplx.__mul(a,b)
	a = cplx.from(a)
	b = cplx.from(b)
	return cplx(
		a[1] * b[1] - a[2] * b[2],
		a[1] * b[2] + a[2] * b[1]
	)
end

function cplx.__div(a,b)
	a = cplx.from(a)
	b = cplx.from(b)
	return a * b:conj() * (1 / b:lenSq())
end

function cplx.__pow(a,b)
	a = cplx.from(a)
	b = cplx.from(b)
	return (a:log() * b):exp()
end

cplx.pow = cplx.__pow

function cplx:conj()
	return cplx(self[1], -self[2])
end

-- aka norm^2, aka modulus^2
function cplx:lenSq()
	return self[1]*self[1] + self[2]*self[2]
end

function cplx:abs()
	return math.sqrt(self:lenSq())
end

function cplx:arg()
	return math.atan2(self[2], self[1])
end

-- returns the first root.  the second is negative the first.
function cplx:sqrt()
	return self^(1/2)
end

function cplx:cbrt()
	return self^(1/3)
end

function cplx:log()
	return cplx(math.log(self:abs()), self:arg())
end

function cplx:exp()
	return math.exp(self[1]) * cplx(math.cos(self[2]), math.sin(self[2]))
end

function cplx:__tostring()
	return self[1]..'+i*'..self[2]
end

cplx.i = cplx(0,1)

return cplx
