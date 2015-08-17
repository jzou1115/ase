package genome;

public class GenomicCoordinate {
	private final int chromosome;
	private final long coord;
	
	public GenomicCoordinate(int i, long newIndex){
		chromosome = i;
		coord = newIndex;
	}
	
	public long distance(GenomicCoordinate other){
		if(chromosome != other.chromosome){
			throw new RuntimeException("Cannot find the distance between coordinates on different chromosomes");
		}
		return Math.abs(this.getCoord() - other.getCoord());
	}
	
	public GenomicCoordinate increment(int numBases){
		long newIndex = coord + numBases;
		return new GenomicCoordinate(chromosome, newIndex);
	}
	
	public GenomicCoordinate decrement(int numBases){
		long newIndex = coord - numBases;
		if(newIndex < 1){
			return new GenomicCoordinate(chromosome, 0);
		}
		return new GenomicCoordinate(chromosome, newIndex);
	}
	
	public int getChromosome(){
		return chromosome;
	}
	
	public long getCoord(){
		return coord;
	}
	
	@Override
	public String toString(){
		return chromosome+"\t"+coord;
	}
	
	@Override
	public boolean equals(Object o){
		if(o == null) return false;
		if(this == o) return true;
		if(o instanceof GenomicCoordinate){
			GenomicCoordinate other = (GenomicCoordinate) o;
			return chromosome == other.chromosome && coord == other.coord;
		}
		return false;
	}

	public int compareTo(GenomicCoordinate o){
		int chrComp = chromosome - o.chromosome;
		if(chrComp != 0) return chrComp;
		if(coord > o.coord){
			return 1;
		}
		if(coord < o.coord){
			return -1;
		}
		return 0;
	}

}
