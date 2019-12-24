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
local cplx = require 'cplx'

-- poly coeffs
local coeffs = table{[0]=-1, 0, 1}:map(function(x) return cplx(x) end)
coeffs.n = #coeffs	-- only exists for the sake of gui input

--[[
local var = symmath.var
local x = var'x'
local f = coeffs[0]
local xi = x
for i=1,coeffs.n do
	f = f + coeffs[i] * xi
	xi = (xi * x)()
end
f = f()
print(var'f':eq(f))
local fFunc = symmath.export.Lua:compile(f, {x})
--]]
local function fFunc(x)
	local sum = coeffs[0]
	local xi = x
	for i=1,coeffs.n do
		sum = sum + xi * coeffs[i]
		xi = xi * x
	end
	return sum
end
		
local sqrt3 = math.sqrt(3)

local function fRootsAt(y)
	local n = coeffs.n
	while n > 0 do
		if coeffs[n] ~= 0 then break end
		n = n - 1
	end

	if n == 1 then
		local b, a = table.unpack(coeffs, 0, 1)
		b = b - y
		return {-b / a}
	elseif n == 2 then
		local c, b, a = table.unpack(coeffs, 0, 2)
		c = c - y
		local sqrtD = cplx.sqrt(b*b - 4*a*c)
		return {
			(-b + sqrtD) / (2 * a),
			(-b - sqrtD) / (2 * a),
		}
	elseif n == 3 then
		local d, c, b, a = table.unpack(coeffs, 0, 3)
		d = d - y
		local a2 = b/a
		local a1 = c/a
		local a0 = d/a
		-- http://mathworld.wolfram.com/CubicFormula.html
		local R = (9 * a1 * a2 - 27 * a0 - 2 * a2^3) / 54
		local Q = (3 * a1 - a2^2) / 9
		local D = Q^3 + R^2
		local sqrtD = cplx.sqrt(D)
		local S = cplx.cbrt(R + sqrtD)
		local T = cplx.cbrt(R - sqrtD)
		local B = S + T
		local A = S - T
		return {
			-a2/3 + B,
			-a2/3 - B/2 + cplx.i * sqrt3/2 * A,
			-a2/3 - B/2 - cplx.i * sqrt3/2 * A,
		}
	elseif n == 4 then
do return {} end		
		local e, d, c, b, a = table.unpack(coeffs, 0, 4)
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

	inputInt('n', coeffs, 'n')
	for i=0,coeffs.n do
		coeffs[i] = coeffs[i] or cplx()
		ig.igText('c'..i)
		ig.igSameLine()
		ig.igPushIDStr('c'..i..'re')
		coeffChanged = inputFloat('', coeffs[i], 1) or coeffChanged
		ig.igSameLine()
		ig.igPopID()
		ig.igPushIDStr('c'..i..'im')
		coeffChanged = inputFloat('', coeffs[i], 2) or coeffChanged
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

	gl.glBegin(gl.GL_LINES)
	gl.glColor3f(.1, .1, .1)
	for x=math.floor(xmin), math.ceil(xmax) do
		gl.glVertex2f(x, ymin)
		gl.glVertex2f(x, ymax)
	end
	for y=math.floor(ymin), math.ceil(ymax) do
		gl.glVertex2f(xmin, y)
		gl.glVertex2f(xmax, y)
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
		self.rootsPts = table()
		for j=1,self.height do
			local v = (j-.5) / self.height
			local y = v * ymax + (1 - v) * ymin
			local roots = fRootsAt(y)
			for _,root in ipairs(roots) do
				root[3] = y
				self.rootsPts:insert(root)
			end
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
		gl.glVertex2f(x, fFunc(x)[1])	-- fFunc is cplx ... I'm just looking at the real value ... though I could just as well only look at the magnitude?
	end
	gl.glEnd()
	--]]
	
	App.super.update(self)
end

App():run()
