local class = require 'ext.class'
local math = require 'ext.math'

local Complex = class()

function Complex:init(re,im)
	if type(re) == 'table' and Complex:isa(re) then
		self[1] = re[1] or 0
		self[2] = re[2] or 0
	else
		self[1] = re or 0
		self[2] = im or 0
	end
end

function Complex.from(a)
	if type(a) == 'table' then
		if Complex:isa(a) then
			return a
		else
			return Complex(a)
		end
	elseif type(a) == 'number' then
		return Complex(a,0)
	end
	return error("couldn't convert type "..type(a).." to Complex")
end

function Complex:clone()
	return Complex(self[1], self[2])
end

function Complex.__unm(a)
	a = Complex.from(a)
	return Complex(-a[1], -a[2])
end

function Complex.__add(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return Complex(a[1] + b[1], a[2] + b[2])
end

function Complex.__sub(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return Complex(a[1] - b[1], a[2] - b[2])
end

function Complex.__mul(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return Complex(
		a[1] * b[1] - a[2] * b[2],
		a[1] * b[2] + a[2] * b[1]
	)
end

function Complex.__div(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return a * b:conj() * (1 / b:lenSq())
end

function Complex.__pow(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return (a:log() * b):exp()
end

Complex.pow = Complex.__pow

function Complex.__eq(a,b)
	a = Complex.from(a)
	b = Complex.from(b)
	return a[1] == b[1] and a[2] == b[2]
end

function Complex:conj()
	return Complex(self[1], -self[2])
end

-- aka norm^2, aka modulus^2
function Complex:lenSq()
	return self[1]*self[1] + self[2]*self[2]
end

function Complex:abs()
	return math.sqrt(self:lenSq())
end

function Complex:arg()
	return math.atan2(self[2], self[1])
end

-- returns the first root.  the second is negative the first.
function Complex:sqrt()
	return self^(1/2)
end

function Complex:cbrt()
	return self^(1/3)
end

function Complex:log()
	return Complex(math.log(self:abs()), self:arg())
end

function Complex:exp()
	return math.exp(self[1]) * Complex(math.cos(self[2]), math.sin(self[2]))
end

function Complex:__tostring()
	if self[2] == 0 then
		return tostring(self[1])
	elseif self[1] == 0 then
		return 'i*'..self[2]
	end
	return self[1]..'+i*'..self[2]
end

function Complex:isfinite()
	return math.isfinite(self[1]) and math.isfinite(self[2])
end

Complex.i = Complex(0,1)

return Complex
