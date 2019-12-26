#!/usr/bin/env luajit

require 'ext'
local symmath = require 'symmath'
local gl = require 'gl'
local ig = require 'ffi.imgui'
local ffi = require 'ffi'
local Orbit = require 'glapp.orbit'
local View = require 'glapp.view'
local ImGuiApp = require 'imguiapp'
local vec3 = require 'vec.vec3'
local quat = require 'vec.quat'
local Complex = require 'complex'
--local Complex = require 'symmath.complex'	-- hmm
local Poly = require 'poly'

-- polynomial function, in terms of coefficients
local poly = Poly{[0]=0, 0, 1}

local sqrt3 = math.sqrt(3)

--local zExpr = var'x' + symmath.i * var'y'
--local fExpr = poly[0] + 
-- [[ root finding attempt based on newton descent over the function norm sq
-- f(x) = c_i x^i
-- |f(x)|^2 = f(x) f*(x) = (c_i x^i) (c_j x^j)* = (c_i x^i) ((c_j)* (x*)^j)
-- df/dx_re = i c_i,re x^(i-1) re
local function buildFindRootFuncs()
	local symmath = require 'symmath'
	local var = symmath.var
	local Matrix = symmath.Matrix
	local x = var'x'
	local y = var'y'
	local yc = var'yc'
	local i = var'i'
	local z = x + i * y
	local zk = 1
	local f = poly[0][1] - yc + i * poly[0][2]
	for k=1,poly.n do
		zk = zk * z
		local ck = poly[k][1] + i * poly[k][2]
		f = f + ck * zk
		f = f():replace(i^2, -1)()
	end
	--[=[ using gradient descent / optionally Hessian inverse on the real parameters of the complex function's modulus squared
	--f * conj(f)
	f = (f * f:replace(i, -i))():replace(i^2, -1)()
	local df_dx = f:diff(x)()
	local df_dy = f:diff(y)()
	local dz_dn = Matrix({df_dx}, {df_dy})
	-- [[
	local d2f_dxx = df_dx:diff(x)()
	local d2f_dyy = df_dy:diff(y)()
	local d2f_dxy = df_dx:diff(y)()
	dz_dn = Matrix({d2f_dxx, d2f_dxy}, {d2f_dxy, d2f_dyy}):inverse() * dz_dn
	--]]
	print(dz_dn)
	print()
	dz_dn = dz_dn()
	print(dz_dn)
	print()
	dz_dn = table{dz_dn[1][1], dz_dn[2][1]}
	--]=]
	-- [=[ using Newton descent
	-- optionally to get away from the real axis: g(z) = exp(i epsilon) f(z), so z = z - f(z)/f'(z) becomes z = z - f(z) / (f'(z) + i epsilon f(z))
	local df_dz = poly[1][1] + i * poly[1][2]
	local zkMinus1 = 1
	for k=2,poly.n do
		zkMinus1 = zkMinus1 * z
		local ck = poly[k][1] + i * poly[k][2]
		df_dz = df_dz + k * ck * zkMinus1
		df_dz = df_dz():replace(i^2, -1)()
	end
	print('f\n'..f)
	print('df/dz\n'..df_dz)
	local dz_dn = -f / df_dz	--newton update: dz/dn
	if symmath.op.div.is(dz_dn) then	-- get rid of i's in the denominator
		local denom = dz_dn[2]:clone()
		dz_dn[1] = (dz_dn[1] * denom:replace(i, -i))()
		dz_dn[2] = (dz_dn[2] * denom:replace(i, -i))()
	end
	dz_dn = dz_dn():replace(i^2, -1)()
	print('dz/dn\n'..dz_dn)
	local poly = dz_dn:polyCoeffs(i)
	local coeffKeys = table.keys(poly):sort()
	assert(#coeffKeys <= 2)
	assert(coeffKeys[1] == 0 or coeffKeys[1] == 1)
	assert(coeffKeys[2] == 0 or coeffKeys[2] == 1)
	dz_dn = {
		poly[0] or 0,	-- coeff of 1
		poly[1] or 0,	-- coeff of i
	}
	
	--]=]

	local dx_dn_func = dz_dn[1]:compile{x,y,yc}
	local dy_dn_func = dz_dn[2]:compile{x,y,yc}
	return function(z, y)
		return Complex(
			dx_dn_func(z[1], z[2], y),
			dy_dn_func(z[1], z[2], y))
	end
end

local function useFindRootFuncs(dz_dt, z, y)
--	print('z', z)
	local maxiter = 1000
	for iter=1,maxiter do
		print('z', z, 'y', y)
		local dz = dz_dt(z, y)
		print('dz', dz)
		if not math.isfinite(z[1]) or not math.isfinite(z[2]) then
			error("got a nan value")
		end
		if dz:lenSq() < 1e-7 then return z end
		z = z - dz
--		print('z', z)
	end
	print("maxiter reached")
end
--]]

local function fRootsAt(y, lastZs)
	local n = poly.n
	while n > 0 do
		if poly[n] ~= Complex(0,0) then break end
		n = n - 1
	end

-- [[ numeric root finding ... using newton's method
	local results = table()
	
	-- but what if all the seed points chosen still converge to the same basin?
	-- the only way to get around this is if you divide out the previously-found roots with polynomial division
	for _,z0 in ipairs(lastZs) do
		for _,z0ofs in ipairs{Complex(1,0), Complex(0,1), Complex(-1,0), Complex(0,-1)} do

			local solvePoly = poly - y
			local solvePolyDiff = solvePoly:diff()
			
			local z0epsilon = 1e-3
			local z = z0 + z0ofs * z0epsilon

			local found
			local maxiter = 100
			for j=1,maxiter do
				local dz_dj = -solvePoly(z) / solvePolyDiff(z)
				if 
				not math.isfinite(dz_dj[1])
				or not math.isfinite(dz_dj[2])
				or dz_dj:lenSq() < 1e-7 
				then 
					found = true
					break
				end
				z = z + dz_dj
			end
			if found then
				if not results:find(nil, function(prevz) return (z - prevz):lenSq() < 1e-7 end) then
					results:insert(z)
				end
				
				-- TODO here ... polynomial long division on (z - z0) from solvePoly = p(z)
			end
		end
	end
	do return results end
--]]
	
	if n == 1 then
		local b, a = table.unpack(poly, 0, 1)
		b = b - y
		return {-b / a}
	elseif n == 2 then
		local c, b, a = table.unpack(poly, 0, 2)
		c = c - y
		local sqrtD = Complex.sqrt(b*b - 4*a*c)
		return {
			(-b + sqrtD) / (2 * a),
			(-b - sqrtD) / (2 * a),
		}
	elseif n == 3 then
		local d, c, b, a = table.unpack(poly, 0, 3)
		d = d - y
		local a2 = b/a
		local a1 = c/a
		local a0 = d/a
		-- http://mathworld.wolfram.com/CubicFormula.html
		local R = (9 * a1 * a2 - 27 * a0 - 2 * a2^3) / 54
		local Q = (3 * a1 - a2^2) / 9
		local D = Q^3 + R^2
		
		local root2s = table{Complex(1,0), Complex(-1,0)}
		local root3s = table{Complex(1,0), Complex(-.5, sqrt3), Complex(-.5, -sqrt3)}

		local results = table()
		do local sqrtDsign = 1 
		--for _,sqrtDsign in ipairs(root2s) do
			local sqrtD = Complex.sqrt(D) * sqrtDsign 
			do local Ssign = 1
			--for _,Ssign in ipairs(root3s) do
				local S = Complex.cbrt(R + sqrtD) * Ssign
				do local Tsign = 1
				--for _,Tsign in ipairs(root3s) do
					local T = Complex.cbrt(R - sqrtD) * Tsign
					local B = S + T
					local A = S - T
					results:insert(-a2/3 + B)
					results:insert(-a2/3 - B/2 + Complex.i * sqrt3/2 * A)
					results:insert(-a2/3 - B/2 - Complex.i * sqrt3/2 * A)
				end
			end
		end
		return results
	elseif n == 4 then
do return {} end		
		local e, d, c, b, a = table.unpack(poly, 0, 4)
		e = e - y

		local a2 = a * a
		local a3 = a * a2
		local b2 = b * b
		local b3 = b * b2
		local b4 = b2 * b2
		local c2 = c * c
		local c3 = c * c2
		local c4 = c2 * c2
		local d2 = d * d
		local d3 = d * d2
		local d4 = d2 * d2
		local e2 = e * e
		local e3 = e * e2
		
		local Delta0 = c2 - 3 * b * d + 12 * a * e
		local Delta1 = 2 * c3 - 9 * b * c * d + 27 * b2 * e + 27 * a * d2 - 72 * a * c * e

--		local Delta = 256 * a3 * e3 - 192 * a2 * b * d * e2 - 128 * a2 * c2 * e2 + 144 * a2 * c * d2 * e - 27 * a2 * d4  
--					+ 144 * a * b2 * c * e2 - 6 * a * b2 * d2 * e - 80 * a * b * c2 * d * e + 18 * a * b * c * d3 + 16 * a * c4 * e 
--					- 4 * a * c3 * d2 - 27 * b4 * e2 + 18 * b3 * c * d * e - 4 * b3 * d3 - 4 * b2 * c3 * e + * b2 * c2 * d2
		
		local Q = math.cbrt((Delta1 + math.sqrt(Delta1^2 - 4 * Delta0^3))/2)
		local p = (8 * a * c - 3 * b2) / (8 * a2)
		local q = (b3 - 4 * a * b * c + 8 * a2 * d) / (8 * a2)
		local S = .5 * math.sqrt(-2/3*p + (Q + Delta0/Q)/(3*a))	
	
		local P = 8 * a * c - 3 * b2
		local R = b3 + 8 * a2 * d - 4 * a * b * c
		local D = 64 * a3 * e - 16 * a2 * c2 + 16 * a * b2 * c - 16 * a2 * b * d - 3 * b4

	end
	return {}
end

local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'roots of polynomials'

function App:init(...)
	App.super.init(self, ...)
	self.view.ortho = true
end

local int = ffi.new('int[1]', 0)
local function inputInt(label, t, k)
	int[0] = tonumber(t[k]) or 0
	if ig.igInputInt(label, int, 1, 10, 0) then
		t[k] = int[0]
		return true
	end
end

local float = ffi.new('float[1]', 0)
local function inputFloat(label, t, k)
	float[0] = tonumber(t[k]) or 0
	if ig.igInputFloat(label, float, .03, .3, '%f', 0) then
		t[k] = float[0]
		return true
	end
end


local coeffChanged
local projecti = false
function App:updateGUI()
	if ig.igButton(self.view.ortho and 'ortho' or 'frustum') then
		self.view.ortho = not self.view.ortho
		if self.view.ortho then
			self.view.angle = quat(0,0,0,1)
			self.view.pos = vec3(0,0,10)	-- = self.view.pos0[3]
		end
	end

	if ig.igButton('project i:'..(projecti and 'yes' or 'no')) then
		projecti = not projecti  
	end

	inputInt('n', poly, 'n')
	for i=0,poly.n do
		poly[i] = poly[i] or Complex()
		ig.igText('c'..i)
		ig.igSameLine()
		ig.igPushIDStr('c'..i..'re')
		coeffChanged = inputFloat('', poly[i], 1) or coeffChanged
		ig.igSameLine()
		ig.igPopID()
		ig.igPushIDStr('c'..i..'im')
		coeffChanged = inputFloat('', poly[i], 2) or coeffChanged
		ig.igPopID()
	end
end

function App:update()
	self.view:setup(self.width/self.height)
	
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)

	local xmin, xmax, ymin, ymax
	if self.view.ortho then
		xmin, xmax, ymin, ymax = self.view:getBounds(self.width / self.height)
		xmin = xmin + self.view.pos[1]
		xmax = xmax + self.view.pos[1]
		ymin = ymin + self.view.pos[2]
		ymax = ymax + self.view.pos[2]
	else
		xmin, xmax, ymin, ymax = -10, 10, -10, 10
	end

	local xsize = xmax - xmin
	local ysize = ymax - ymin

	-- how to handle frustum display / complex axis? keep it proportional to the real axis, but centered at zero
	local zsize = xsize
	local zmax = zsize/2
	local zmin = -zmax

	local mins = vec3(xmin, ymin, zmin)
	local maxs = vec3(xmax, ymax, zmax)

	gl.glBegin(gl.GL_LINES)
	gl.glColor3f(.1, .1, .1)
	for i=1,3 do
		local j = i%3+1
		local k = j%3+1
		for xi=math.floor(mins[i]), math.ceil(maxs[i]) do
			local v = vec3()
			v[i] = xi
			v[j] = mins[j]
			gl.glVertex3f(v[1], v[2], v[3])
			v[j] = maxs[j]
			gl.glVertex3f(v[1], v[2], v[3])
		end
		for xj=math.floor(mins[j]), math.ceil(maxs[j]) do
			local v = vec3()
			v[j] = xj
			v[i] = mins[i]
			gl.glVertex3f(v[1], v[2], v[3])
			v[i] = maxs[i]
			gl.glVertex3f(v[1], v[2], v[3])
		end
	end

	gl.glColor3f(.5, .5, .5)
	gl.glVertex2f(xmin, 0)
	gl.glVertex2f(xmax, 0)
	gl.glVertex2f(0, ymin)
	gl.glVertex2f(0, ymax)
	gl.glEnd()

	if ymin ~= self.cachedYMin
	or ymax ~= self.cachedYMax
	or coeffChanged
	then
		self.cachedYMin = ymin
		self.cachedYMax = ymax
		coeffChanged = false

		--[[ numeric solution
		local dz_dt = buildFindRootFuncs()
		local z0 = Complex(math.random()*2-1, math.random()*2-1)
		--]]

		self.rootsPts = table()
		local lastRoots = {Complex()}
		for j=1,self.height do
			local v = (j-.5) / self.height
			local y = v * ymax + (1 - v) * ymin
		
			--[[ numeric solution
			local root = useFindRootFuncs(dz_dt, z0, y)
--			print(y, root)
			self.rootsPts:insert(root)
			z0 = root	-- use the last root as the seed for the next root
			--]]	
			-- [[ exact solution
			local roots = fRootsAt(y, lastRoots)
			for _,root in ipairs(roots) do
				root[3] = y
				self.rootsPts:insert(root)
			end
			lastRoots = roots
			--]]	
		end
	end

	-- [[ here is a graph of [roots of f(x)-c, c]
	-- I'll have to do it as points for now since there could be multiple roots, and I'm not tracking information on which root of one y is related to which root of another y
	gl.glColor3f(1,1,0)
	gl.glPointSize(3)
	gl.glBegin(gl.GL_POINTS)
	for _,root in ipairs(self.rootsPts) do
		local y = root[3]
		if self.view.ortho then
			gl.glVertex2f(root[1] + root[2] * (projecti and 1 or 0), y)
		else
			gl.glVertex3f(root[1], y, root[2])
		end
	end
	gl.glEnd()
	gl.glPointSize(1)
	--]]

	-- [[ here is a graph of the function: [x,f(x)]
	gl.glColor3f(0,1,0)
	gl.glBegin(gl.GL_LINE_STRIP)
	for i=1,self.width do
		local u = (i-.5) / self.width
		local x = u * xmax + (1 - u) * xmin
		gl.glVertex2f(x, poly(x)[1])	-- poly is complex ... I'm just looking at the real value ... though I could just as well only look at the magnitude?
	end
	gl.glEnd()
	--]]
	
	App.super.update(self)
end

App():run()
