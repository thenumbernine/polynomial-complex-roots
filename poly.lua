local class = require 'ext.class'
local table = require 'ext.table'
local Complex = require 'complex'

local Poly = class()

-- poly is indexed 0..self.n
function Poly:init(t)
	self.n = table.maxn(t)
	for i=0,self.n do
		self[i] = Complex(t[i] or 0)
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

return Poly
