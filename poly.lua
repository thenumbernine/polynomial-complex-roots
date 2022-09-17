local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local Complex = require 'complex'

local Poly = class()

-- poly is indexed 0..self.n
function Poly:init(t)
	if t == nil then
		self.n = 0
		self[0] = Complex()
	elseif type(t) == 'number' then
		self.n = 0
		self[0] = Complex.from(t)
	elseif type(t) == 'table' then
		if Complex:isa(t) then
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

function Poly.from(x)
	if Poly:isa(x) then return x end
	if type(x) == 'number' or Complex:isa(x) then
		return Poly{[0]=x}
	end
	return Poly(x)
end

function Poly.fromRoot(z)
	return Poly{[0] = -z, 1}
end

function Poly:__call(z)
	local sum = self[0] or Complex()
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
	for k=self.n,0,-1 do
		if self[k] ~= Complex(0,0) then
			local coeffstr = ''
			if k == 0 then
				coeffstr = tostring(self[k])
				if coeffstr:find'%+' then coeffstr = '('..coeffstr..')' end
				if result ~= '' then result = result .. ' + ' end
				result = result .. coeffstr
			else
				if self[k] ~= Complex(1,0) then
					coeffstr = tostring(self[k])
					if coeffstr:find'%+' then coeffstr = '('..coeffstr..')' end
				end
				if k == 1 then
					if coeffstr ~= '' then coeffstr = coeffstr..' * ' end
					if result ~= '' then result = result .. ' + ' end
					result = result .. coeffstr..'z'
				else
					if coeffstr ~= '' then coeffstr = coeffstr..' * ' end
					if result ~= '' then result = result .. ' + ' end
					result = result .. coeffstr..'z^'..k
				end
			end
		end
	end
	if result == '' then return '0' end
	return result
end

function Poly:removeExtraZeroes()
	while self.n > 0 and self[self.n] == Complex(0,0) do
		self[self.n] = nil
		self.n = self.n - 1
	end
end

function Poly.__add(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	local c = Poly()
	c.n = math.max(a.n, b.n)
	for k=0,c.n do
		c[k] = (a[k] or 0) + (b[k] or 0)
	end
	c:removeExtraZeroes()
	return c
end

function Poly.__sub(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	local c = Poly()
	c.n = math.max(a.n, b.n)
	for k=0,c.n do
		c[k] = (a[k] or Complex()) - (b[k] or Complex())
	end
	c:removeExtraZeroes()
	return c
end

function Poly.__mul(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	local c = Poly()
	c.n = a.n + b.n
	for i=0,c.n do
		-- c_i = sum j,k such that j+k == i of a_j*b_k
		local sum = 0
		for j=0,math.min(i,a.n) do
			local k = i - j
			if k <= b.n then
				sum = sum + a[j] * b[k]
			end
		end
		c[i] = sum
	end
	return c
end

function Poly:lshift(amount)
	local result = Poly()
	result.n = self.n + amount
	for i=0,self.n do
		result[i+amount] = Complex(self[i])
	end
	for i=0,amount-1 do
		result[i] = Complex()
	end
	return result
end

function Poly.div(a,b)
	a = Poly.from(a)
	b = Poly.from(b)
	
	--[[ polynomial long division
	a / b
	= (a_m z^n + ... a_1 z + a_0) / (b_n z^n + ... b_1 z + b_0)
	
	                               ____________________________
	= (b_n z^n + ... b_1 z + b_0) | (a_m z^n + ... a_1 z + a_0)
	--]]

	local c = Poly()
	for k=a.n,b.n,-1 do
--print('--------------------------')
--print('k = '..k)		
--print('k - b.n = '..(k - b.n))
--print('a[k] = '..a[k])
--print('b[b.n] = '..b[b.n])
		if a[k] then	-- otherwise a[k] == 0 and we don't need this step (but don't forget to scale up the denominator)
			local coeff = a[k] / b[b.n]
--print('coeff = '..coeff)
			if coeff ~= Complex(0,0) then
--print('a = '..a)
--print('b = '..b)
--print('b * coeff = '..(b * coeff))
--print('(b * coeff):lshift(k - b.n) = '..((b * coeff):lshift(k - b.n)))
--print('a - (b * coeff):lshift(k - b.n) = '..(a - (b * coeff):lshift(k - b.n)))
				a = a - (b * coeff):lshift(k - b.n)
--print('next a = '..a)
				c[k - b.n] = coeff
--print('setting c['..(k - b.n)..'] = '..c[k - b.n])
				c.n = math.max(c.n, k - b.n)
				a[k] = Complex()
			end
		else
			c[k - b.n] = Complex()
		end
	end
	a:removeExtraZeroes()
--print('--------------------------')
--print('result = '..c)
--print('remainder = '..a)
	return c, a
end

Poly.__div = Poly.div

--[[ test code
local a = Poly{[0]=-1, 1, 1, 1}
print('a = '..a)
local b = Poly{[0]=-.54369, 1}
--local b = Poly{[0]=-1, 1}
print('b = '..b)
print('a*b = '..(a*b))
print('a/b = '..table{a:div(b)}:mapi(function(x) return tostring(x) end):concat' R ')
os.exit()
--]]

return Poly
