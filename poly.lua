local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local Complex = require 'complex'

local Poly = class()

-- poly is indexed 0..self.n
function Poly:init(t)
	if not t then
		self.n = 0
	else
		self.n = table.maxn(t)
		for k=0,self.n do
			self[k] = Complex(t[k] or 0)
		end
	end
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
				coeffstr = coeffstr .. ' * '
			end
			if k == 0 then
				result = result .. sep .. coeffstr
			elseif k == 1 then
				result = result .. sep .. coeffstr..'z'
			else
				result = result .. sep .. coeffstr..'z^'..k
			end
			sep = ' + '
		end
	end
	if result == '' then return '0' end
	return result
end

return Poly
