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
