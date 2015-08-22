#!/usr/bin/ruby

require 'pp'
include Math

class CubicLattice
  def initialize max
    @max = max
    @coord = [-1,0,0]
  end
  def increment i
    @coord[i] += 1
    if @coord[i] >= @max
      @coord[i] = 0
      increment i+1
    end
  end
  def next
    increment 0
    @coord
  end
end

class Simulator
  attr_accessor :box_size

  def initialize(energy_cal, eq_move, pro_move)
    @equilibrate_move = eq_move
    @production_move = pro_move

    @n_atoms = 100  # number of atoms
    @box_size = 40  # box size (amstrong)

    @energy_cal = energy_cal.new(self)

    @max_translate = 0.4  # atom max translation angstroms

    @kT = 1.987206504191549e-3 * 298.15   # Simulation temperature * k (kcal mol-1 K-1)

    @n_accept = 0
    @n_reject = 0

    # put particles in lattice
    @particles = []

    atom_per_axis = (@n_atoms**(1.0/3)).ceil
    spacing = @box_size.to_f / atom_per_axis

    coord = CubicLattice.new(atom_per_axis)
    @n_atoms.times do
      @particles.push coord.next.map {|c| 0.01 + c*spacing }
    end

    # singleton methods for enumerating particles
    def @particles.each_pair
      (0...length).each do |i|
        ((i+1)...length).each do |j|
          yield [i,j]
        end
      end
    end

    def @particles.each_pair_with_i(iatom)
      (0...length).each {|j| yield [iatom, j].sort if iatom != j }
    end

    @energy_grid = {}
    @energy_cal.updateall(@energy_grid, @particles)
    @energy = @energy_grid.each_value.reduce :+

  end

  def run
    puts "perform #@equilibrate_move equilibration"
    @equilibrate_move.times {try_move }

    puts "perform #@production_move steps for production"
    print "%-13s: %-10s %-10s\n" % ['move number', 'energy', 'acceptance']

    accumulate_vol = 0
    @production_move.times do |n_trial|
      try_move
      accumulate_vol += @box_size**3
      if n_trial % 10 == 0
        print "%-13d: %-10.2f %-10s\n" % [n_trial, @energy, (@n_accept*100/(@n_reject+@n_accept)).to_s + '%']
      end
      # print_pdb n_trial if n_trial % 100 == 0
    end
    puts "=" * 40
    puts "mean volume: #{accumulate_vol/@production_move}"
  end

  def try_move
    if test_move
      @n_accept += 1
    else
      @n_reject += 1
      restore_state
    end
  end

  def test_move
    raise "test_move() has not overriden by subclass"
  end

  def translation_move
    i = rand(@n_atoms)
    @backup = [:trans, i, atom_backup(i)]
    @particles[i].each_index do |j|
      @particles[i][j] += rand((-@max_translate)..@max_translate)
      if not (0..@box_size).cover? @particles[i][j]
        @particles[i][j] %= @box_size
      end
    end
    @energy_cal.update_ith(@energy_grid, @particles, i)
    @energy = @energy_grid.each_value.reduce :+
  end

  def atom_backup iatom
    e = {}
    @particles.each_pair_with_i(iatom) {|pair| e[pair] = @energy_grid[pair] }
    [@particles[iatom].dup, e]
  end

  def restore_state
    case @backup[0]
    when :trans
      iatom, state = @backup[1], @backup[2]
      @particles[iatom] = state[0]
      state[1].each_pair{|key, value| @energy_grid[key] = value }
    when :vol
      scale = @backup[1]
      @particles.each do |p|
        p.each_index {|i| p[i] /= scale }
      end
      @box_size /= scale
    end
    @energy_cal.updateall(@energy_grid, @particles)
    @energy = @energy_grid.each_value.reduce :+
  end

  def print_pdb(nth)
    open("output%06d.pdb" % nth, 'w') do |out|
      out.printf("CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n", @box_size, @box_size, @box_size)
      @particles.each_with_index do |p, i|
        out.printf("ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                   i, p[0], p[1], p[2])
        out.puts "TER"
      end
    end
  end
end

class SimulatorNVT < Simulator
  def test_move
    old_E = @energy
    translation_move
    new_E = @energy
    new_E <= old_E or exp(-(1/@kT)*(new_E - old_E)) >= rand
  end

end

class SimulatorNPT < Simulator
  attr_accessor :pressure

  def test_move
    old_E = @energy
    type = rand < 0.9 ? :trans : :vol
    if type == :trans
      translation_move
    else
      old_vol = @box_size ** 3
      volume_move
      new_vol = @box_size ** 3
    end
    new_E = @energy

    vol_factor = (type == :trans) ? 0 : -(1/@kT)*@pressure*(new_vol-old_vol) + (@n_atoms+1)*log(new_vol/old_vol)
    # puts "new_E: #{new_E} old_E: #{old_E} vol_factor: #{vol_factor}"

    tmp = exp(-(1/@kT)*(new_E - old_E) + vol_factor) >= rand
    if type != :trans
      if tmp == true
        puts "accept volume move: deltaE #{new_E-old_E}, vol_factor: #{vol_factor}"
      else
        # puts "exp(-(1/@kT)*(new_E - old_E) + vol_factor) = #{exp(-(1/@kT)*(new_E - old_E) + vol_factor)}"
        # puts "(@n_atoms+1)*log(new_vol/old_vol) = #{(@n_atoms+1)*log(new_vol/old_vol)}"
        # puts "-(1/@kT)*@pressure*(new_vol-old_vol) = #{-(1/@kT)*@pressure*(new_vol-old_vol)}"
        puts "reject volume move: deltaE #{new_E-old_E}, vol_factor: #{vol_factor}"
      end
    else
      if tmp == true
        puts "accept translation move"
      else
        puts "reject translation move"
      end
    end
    tmp
  end

  def volume_move
    delta_vol = 1
    @old_size = @box_size
    new_vol = exp(log(@box_size**3) + (rand-0.5) * delta_vol)
    @box_size = new_vol**(1.0/3)
    scale = (@box_size/@old_size)
    puts "scaling: #{scale}"
    @backup = [:vol, scale]

    @particles.each do |p|
      p.each_index {|i| p[i] *= scale }
    end
    @energy_cal.updateall(@energy_grid, @particles)
    @energy = @energy_grid.each_value.reduce :+
  end
end

class MinImageEnergy

  # Give the Lennard Jones parameters for the atoms
  # (these are the OPLS parameters for Krypton)
  Sigma = 3.624     # angstroms
  Epsilon = 0.317   # kcal mol-1

  def initialize host
    @simulator = host
  end

  # central cell particle p1 p2, calculate energy between them
  # only minimum image considered
  def interaction(p1, p2)
    r_sq = p1.zip(p2).reduce(0) {|sum, coords| sum + (min_img(coords[1] - coords[0]))**2 }
    r = sqrt r_sq
    4.0 * Epsilon * ( (Sigma/r)**12 - (Sigma/r)**6 )
  end

  def min_img(distance)
    box_size = @simulator.box_size
    distance = distance.abs
    distance > (0.5*box_size) ? box_size - distance : distance
  end

  def updateall(energy_grid, particles)
    particles.each_pair do |key|
      # puts "calculate interaction of points: #{key[0]} #{key[1]}"
      energy_grid[key] = interaction(particles[key[0]], particles[key[1]])
    end
  end

  def update_ith(energy_grid, particles, iatom)
    raise "invalid atom" if ! iatom.is_a? Integer
    particles.each_pair_with_i iatom do |key|
      # puts "calculate interaction of points: #{key[0]} #{key[1]}"
      energy_grid[key] = interaction(particles[key[0]], particles[key[1]])
    end
  end
end

equilibrate_move = ARGV[0] || '500'
production_move = ARGV[1] || '200'
enssemble = ARGV[2] || 'npt'

simulator_table = {'nvt' => SimulatorNVT, 'npt' => SimulatorNPT}

sim = simulator_table[enssemble].new(MinImageEnergy, equilibrate_move.to_i, production_move.to_i)
if enssemble == 'npt'
  sim.pressure = 1 * 1.458e-5  # kcal mol^-1 A^-3
end

sim.run
