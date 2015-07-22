#!/usr/bin/env ruby

# N is window size, n is count of each base in the window, K is
# alphabet size

# def cwf
#   val = @counts.map { |_, count| fact(count) }.reduce(1, :*)
#   # p @counts.map { |_, count| fact(count) }
#   assert val != 0
#   (1.0 / @window) * log_k(fact(@window) / val)
# end

# def ce
#   val = @counts.map do |_, count|
#     quo = count / @window
#     quo * log_k(quo)
#   end

#   -val.reduce(:+)
# end

Signal.trap("PIPE", "EXIT")

methods = File.join(File.expand_path('~'), 'lib', 'ruby', 'ryan.rb')
require_relative methods

Ryan.req *%w[parse_fasta ffi]

module CMath
  extend FFI::Library
  ffi_lib "#{File.dirname(__FILE__)}/lib/complexity.so"

  attach_function :cwf, [:pointer, :long], :double
  attach_function :ce, [:pointer, :long], :double
end

class Complexity
  attr_accessor :window, :counts

  ALPHABET_SIZE = 4
  OFFSET = 0

  # counts is an array of base counts
  def initialize window, counts
    abort "Window must be > 0" unless window > 0
    @window = window.to_f

    unless counts.all? { |count| count > 0 }
      abort "Counts must all be positive"
    end
    @counts = counts
  end

  def cwf
    mp = FFI::MemoryPointer.new(:long, ALPHABET_SIZE)
    mp.put_array_of_long OFFSET, @counts
    CMath.cwf mp, @window
  end


  def ce
    mp = FFI::MemoryPointer.new(:long, ALPHABET_SIZE)
    mp.put_array_of_long OFFSET, @counts
    CMath.ce mp, @window
  end
end

def full_len_stats infile, outfile
  outfile.puts %w[name start stop cwf ce].join "\t"

  FastaFile.open(infile).each_record do |head, seq|
    len = seq.length
    window = seq.length
    counts = seq.base_counts.values
    comp = Complexity.new window, counts
    start = 1
    stop = len

    outfile.puts [head, start, stop, comp.cwf, comp.ce].join "\t"
  end
end

def window_stats infile, outfile, window, step
  outfile.puts %w[name start stop cwf ce].join "\t"

  FastaFile.open(infile).each_record do |head, seq|
    len = seq.length
    abort "window > len" unless window <= len

    (0 .. len - window).step(step) do |n|
      start = n + 1
      stop = n + window
      this_seq = seq[n .. n + window - 1]
      counts = this_seq.base_counts.values
      comp = Complexity.new window, counts

      outfile.puts [head, start, stop, comp.cwf, comp.ce].join "\t"
    end
  end
end

opts = Trollop.options do
  banner <<-EOS

  Do the thing.

  Options:
  EOS

  opt(:infile, "Input file", type: :string)
  opt(:window, "Size of the window", type: :int,
      default: 150)
  opt(:step, "Number of bases to slide the window each time",
      type: :int, default: 1)
  opt(:full, "Do the full length sequence")
  opt(:outdir, "Output directory", type: :string, default: ".")
end

infile = Ryan.check_file(opts[:infile], :infile)
Ryan.try_mkdir(opts[:outdir])

Trollop::die :window unless 1 < opts[:window]
Trollop::die :step unless 0 < opts[:step]

window = opts[:window]
step = opts[:step]

if opts[:full]
  outf = File.join(opts[:outdir],
                   "#{infile[:base]}.complexity_full_length.txt")

  File.open(outf, "w") do |outfile|
    full_len_stats opts[:infile], outfile
  end
else
  outf = File.join(opts[:outdir],
                   "#{infile[:base]}.complexity_window_#{window}_" +
                   "step_#{step}.txt")

  File.open(outf, "w") do |outfile|
    window_stats opts[:infile], outfile, window, step
  end
end
