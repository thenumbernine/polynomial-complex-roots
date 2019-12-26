local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local Complex = require 'complex'

local Poly = class()

-- poly is indexed 0..self.n
function Poly:init(t)
	if t == nil then
		self.n = 0
	elseif type(t) == 'number' then
		self.n = 0
		self[0] = Complex.from(t)
	elseif type(t) == 'table' then
		if Complex.is(t) then
			self.n = 0
			self[0] = t:clone()
		else
			self.n = table.maxn(t)
			for k=0,self.n do
				self[k] = Complex(t[k] or 0)
			end
		end
	end
end

-- hmm, make this a part of 'class', as the cast operator
-- so 'is' is instanceof, 'from' is cast, and 'init' is ctor
function Poly.from(x)
	if Poly.is(x) then return x end
	return Poly(x)
end

function Poly:__call(z)
	local sum = self[0]
	local zk
	for k=1,self.n do
		zk = zk and zk * z or z
		sum = sum + zk * self[k]
	end
	return sum
end

function Poly:diff()
	local diff = Poly()
	if self.n > 0 then
		for k=1,self.n do
			diff[k-1] = Complex(self[k]) * k
		end
		diff.n = self.n-1
	end
	return diff
end

function Poly:clone()
	return Poly(self)
end

function Poly:__tostring()
	local result = ''
	local sep = ''
	for k=0,self.n do
		if self[k] ~= Complex(0,0) then
			local coeffstr = ''
			if self[k] ~= Complex(1,0) then
				coeffstr = tostring(self[k])
				if coeffstr:find'%+' then coeffstr = '('..coeffstr..')' end
				coeffstr = coeffstr
			end
			if k == 0 then
				result = result .. sep .. coeffstr
			elseif k == 1 then
				result = result .. sep .. coeffstr..' * z'
			else
				result = result .. sep .. coeffstr..' * z^'..k
			end
			sep = ' + '
		end
	end
	if result == '' then return '0' end
	return result
end

function Poly.__add(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	local c = Poly()
	c.n = math.max(a.n, b.n)
	for k=0,c.n do
		c[k] = (a[k] or 0) + (b[k] or 0)
	end
	return c
end

function Poly.__sub(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	local c = Poly()
	c.n = math.max(a.n, b.n)
	for k=0,c.n do
		c[k] = (a[k] or 0) - (b[k] or 0)
	end
	return c
end

return Poly
