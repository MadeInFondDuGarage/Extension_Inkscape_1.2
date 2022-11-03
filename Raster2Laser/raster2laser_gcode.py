'''
# ----------------------------------------------------------------------------
#Ensemble repris et modifier pour la version inkscape 1.23 par MadeInFondDuGarage
#Merci a 305 engenering pour leur extension de de base qui m'a permis de refaie quelque chose de fonctionnel
#Version 1.2_10
# ----------------------------------------------------------------------------
'''

import os
import re


import subprocess
import inkex

#pourla partie png
import collections
import io   # For io.BytesIO
import itertools
import math
import struct
import warnings
import zlib
from array import array
from random import randint

#######Début de la partie lier a la lecture du png
# The PNG signature.
# http://www.w3.org/TR/PNG/#5PNG-file-signature
signature = struct.pack('8B', 137, 80, 78, 71, 13, 10, 26, 10)

# The xstart, ystart, xstep, ystep for the Adam7 interlace passes.
adam7 = ((0, 0, 8, 8),
         (4, 0, 8, 8),
         (0, 4, 4, 8),
         (2, 0, 4, 4),
         (0, 2, 2, 4),
         (1, 0, 2, 2),
         (0, 1, 1, 2))

def adam7_generate(width, height):
    """
    Generate the coordinates for the reduced scanlines
    of an Adam7 interlaced image
    of size `width` by `height` pixels.

    Yields a generator for each pass,
    and each pass generator yields a series of (x, y, xstep) triples,
    each one identifying a reduced scanline consisting of
    pixels starting at (x, y) and taking every xstep pixel to the right.
    """

    for xstart, ystart, xstep, ystep in adam7:
        if xstart >= width:
            continue
        yield ((xstart, y, xstep) for y in range(ystart, height, ystep))


# Models the 'pHYs' chunk (used by the Reader)
Resolution = collections.namedtuple('_Resolution', 'x y unit_is_meter')


def group(s, n):
    return list(zip(* [iter(s)] * n))


def isarray(x):
    return isinstance(x, array)


class Error(Exception):
    def __str__(self):
        return self.__class__.__name__ + ': ' + ' '.join(self.args)


class FormatError(Error):
    """
    Problem with input file format.
    In other words, PNG file does not conform to
    the specification in some way and is invalid.
    """


class ProtocolError(Error):
    """
    Problem with the way the programming interface has been used,
    or the data presented to it.
    """


class ChunkError(FormatError):
    pass


class Reader:
    """
    Pure Python PNG decoder in pure Python.
    """

    def __init__(self, _guess=None, filename=None, file=None):
        """
        The constructor expects exactly one keyword argument.
        If you supply a positional argument instead,
        it will guess the input type.
        Choose from the following keyword arguments:

        filename
          Name of input file (a PNG file).
        file
          A file-like object (object with a read() method).
        bytes
          ``bytes`` or ``bytearray`` with PNG data.

        """
        keywords_supplied = (
            (_guess is not None) +
            (filename is not None) +
            (file is not None))
        if keywords_supplied != 1:
            raise TypeError("Reader() takes exactly 1 argument")

        # Will be the first 8 bytes, later on.  See validate_signature.
        self.signature = None
        self.transparent = None
        # A pair of (len,type) if a chunk has been read but its data and
        # checksum have not (in other words the file position is just
        # past the 4 bytes that specify the chunk type).
        # See preamble method for how this is used.
        self.atchunk = None

        if _guess is not None:
            if isinstance(_guess, str):
                filename = _guess
            elif hasattr(_guess, 'read'):
                file = _guess

        if filename is not None:
            self.file = open(filename, "rb")
        else:
            raise ProtocolError("expecting filename, file or bytes array")

    def chunk(self, lenient=False):
        """
        Read the next PNG chunk from the input file;
        returns a (*type*, *data*) tuple.
        *type* is the chunk's type as a byte string
        (all PNG chunk types are 4 bytes long).
        *data* is the chunk's data content, as a byte string.

        If the optional `lenient` argument evaluates to `True`,
        checksum failures will raise warnings rather than exceptions.
        """

        self.validate_signature()

        # http://www.w3.org/TR/PNG/#5Chunk-layout
        if not self.atchunk:
            self.atchunk = self._chunk_len_type()
        if not self.atchunk:
            raise ChunkError("No more chunks.")
        length, type = self.atchunk
        self.atchunk = None

        data = self.file.read(length)
        if len(data) != length:
            raise ChunkError(
                'Chunk %s too short for required %i octets.'
                % (type, length))
        checksum = self.file.read(4)
        if len(checksum) != 4:
            raise ChunkError('Chunk %s too short for checksum.' % type)
        verify = zlib.crc32(type)
        verify = zlib.crc32(data, verify)
        verify = struct.pack('!I', verify)
        if checksum != verify:
            (a, ) = struct.unpack('!I', checksum)
            (b, ) = struct.unpack('!I', verify)
            message = ("Checksum error in %s chunk: 0x%08X != 0x%08X."
                       % (type.decode('ascii'), a, b))
            if lenient:
                warnings.warn(message, RuntimeWarning)
            else:
                raise ChunkError(message)
        return type, data

    def chunks(self):
        """Return an iterator that will yield each chunk as a
        (*chunktype*, *content*) pair.
        """

        while True:
            t, v = self.chunk()
            yield t, v
            if t == b'IEND':
                break

    def undo_filter(self, filter_type, scanline, previous):
        """
        Undo the filter for a scanline.
        `scanline` is a sequence of bytes that
        does not include the initial filter type byte.
        `previous` is decoded previous scanline
        (for straightlaced images this is the previous pixel row,
        but for interlaced images, it is
        the previous scanline in the reduced image,
        which in general is not the previous pixel row in the final image).
        When there is no previous scanline
        (the first row of a straightlaced image,
        or the first row in one of the passes in an interlaced image),
        then this argument should be ``None``.

        The scanline will have the effects of filtering removed;
        the result will be returned as a fresh sequence of bytes.
        """

        # :todo: Would it be better to update scanline in place?
        result = scanline

        if filter_type == 0:
            return result

        if filter_type not in (1, 2, 3, 4):
            raise FormatError(
                'Invalid PNG Filter Type.  '
                'See http://www.w3.org/TR/2003/REC-PNG-20031110/#9Filters .')

        fu = max(1, self.psize)

        if not previous:
            previous = bytearray([0] * len(scanline))

        fn = (None,
              undo_filter_sub,
              undo_filter_up,
              undo_filter_average,
              undo_filter_paeth)[filter_type]
        fn(fu, scanline, previous, result)
        return result

    def _iter_bytes_to_values(self, byte_rows):
        """
        Iterator that yields each scanline;
        each scanline being a sequence of values.
        `byte_rows` should be an iterator that yields
        the bytes of each row in turn.
        """

        for row in byte_rows:
            yield self._bytes_to_values(row)

    def _bytes_to_values(self, bs, width=None):
        """Convert a packed row of bytes into a row of values.
        Result will be a freshly allocated object,
        not shared with the argument.
        """

        if self.bitdepth == 8:
            return bytearray(bs)
        if self.bitdepth == 16:
            return array('H',
                         struct.unpack('!%dH' % (len(bs) // 2), bs))

        assert self.bitdepth < 8
        if width is None:
            width = self.width
        # Samples per byte
        spb = 8 // self.bitdepth
        out = bytearray()
        mask = 2**self.bitdepth - 1
        shifts = [self.bitdepth * i
                  for i in reversed(list(range(spb)))]
        for o in bs:
            out.extend([mask & (o >> i) for i in shifts])
        return out[:width]

    def _iter_straight_packed(self, byte_blocks):
        """Iterator that undoes the effect of filtering;
        yields each row as a sequence of packed bytes.
        Assumes input is straightlaced.
        `byte_blocks` should be an iterable that yields the raw bytes
        in blocks of arbitrary size.
        """

        # length of row, in bytes
        rb = self.row_bytes
        a = bytearray()
        # The previous (reconstructed) scanline.
        # None indicates first line of image.
        recon = None
        for some_bytes in byte_blocks:
            a.extend(some_bytes)
            while len(a) >= rb + 1:
                filter_type = a[0]
                scanline = a[1: rb + 1]
                del a[: rb + 1]
                recon = self.undo_filter(filter_type, scanline, recon)
                yield recon
        if len(a) != 0:
            # :file:format We get here with a file format error:
            # when the available bytes (after decompressing) do not
            # pack into exact rows.
            raise FormatError('Wrong size for decompressed IDAT chunk.')
        assert len(a) == 0

    def validate_signature(self):
        """
        If signature (header) has not been read then read and
        validate it; otherwise do nothing.
        No signature (empty read()) will raise EOFError;
        An invalid signature will raise FormatError.
        EOFError is raised to make possible the case where
        a program can read multiple PNG files from the same stream.
        The end of the stream can be distinguished from non-PNG files
        or corrupted PNG files.
        """

        if self.signature:
            return
        self.signature = self.file.read(8)
        if len(self.signature) == 0:
            raise EOFError("End of PNG stream.")
        if self.signature != signature:
            raise FormatError("PNG file has invalid signature.")

    def preamble(self, lenient=False):
        """
        Extract the image metadata by reading
        the initial part of the PNG file up to
        the start of the ``IDAT`` chunk.
        All the chunks that precede the ``IDAT`` chunk are
        read and either processed for metadata or discarded.

        If the optional `lenient` argument evaluates to `True`,
        checksum failures will raise warnings rather than exceptions.
        """

        self.validate_signature()

        while True:
            if not self.atchunk:
                self.atchunk = self._chunk_len_type()
                if self.atchunk is None:
                    raise FormatError('This PNG file has no IDAT chunks.')
            if self.atchunk[1] == b'IDAT':
                return
            self.process_chunk(lenient=lenient)

    def _chunk_len_type(self):
        """
        Reads just enough of the input to
        determine the next chunk's length and type;
        return a (*length*, *type*) pair where *type* is a byte sequence.
        If there are no more chunks, ``None`` is returned.
        """

        x = self.file.read(8)
        if not x:
            return None
        if len(x) != 8:
            raise FormatError(
                'End of file whilst reading chunk length and type.')
        length, type = struct.unpack('!I4s', x)
        if length > 2 ** 31 - 1:
            raise FormatError('Chunk %s is too large: %d.' % (type, length))
        # Check that all bytes are in valid ASCII range.
        # https://www.w3.org/TR/2003/REC-PNG-20031110/#5Chunk-layout
        type_bytes = set(bytearray(type))
        if not(type_bytes <= set(range(65, 91)) | set(range(97, 123))):
            raise FormatError(
                'Chunk %r has invalid Chunk Type.'
                % list(type))
        return length, type

    def process_chunk(self, lenient=False):
        """
        Process the next chunk and its data.
        This only processes the following chunk types:
        ``IHDR``, ``PLTE``, ``bKGD``, ``tRNS``, ``gAMA``, ``sBIT``, ``pHYs``.
        All other chunk types are ignored.

        If the optional `lenient` argument evaluates to `True`,
        checksum failures will raise warnings rather than exceptions.
        """

        type, data = self.chunk(lenient=lenient)
        method = '_process_' + type.decode('ascii')
        m = getattr(self, method, None)
        if m:
            m(data)

    def _process_IHDR(self, data):
        # http://www.w3.org/TR/PNG/#11IHDR
        if len(data) != 13:
            raise FormatError('IHDR chunk has incorrect length.')
        (self.width, self.height, self.bitdepth, self.color_type,
          self.compression, self.filter,
          self.interlace) = struct.unpack("!2I5B", data)

        check_bitdepth_colortype(self.bitdepth, self.color_type)

        if self.compression != 0:
            raise FormatError(
                "Unknown compression method %d" % self.compression)
        if self.filter != 0:
            raise FormatError(
                "Unknown filter method %d,"
                " see http://www.w3.org/TR/2003/REC-PNG-20031110/#9Filters ."
                % self.filter)
        if self.interlace not in (0, 1):
            raise FormatError(
                "Unknown interlace method %d, see "
                "http://www.w3.org/TR/2003/REC-PNG-20031110/#8InterlaceMethods"
                " ."
                % self.interlace)

        # Derived values
        # http://www.w3.org/TR/PNG/#6Colour-values
        colormap = bool(self.color_type & 1)
        greyscale = not(self.color_type & 2)
        alpha = bool(self.color_type & 4)
        color_planes = (3, 1)[greyscale or colormap]
        planes = color_planes + alpha

        self.colormap = colormap
        self.greyscale = greyscale
        self.alpha = alpha
        self.color_planes = color_planes
        self.planes = planes
        self.psize = float(self.bitdepth) / float(8) * planes
        if int(self.psize) == self.psize:
            self.psize = int(self.psize)
        self.row_bytes = int(math.ceil(self.width * self.psize))
        # Stores PLTE chunk if present, and is used to check
        # chunk ordering constraints.
        self.plte = None
        # Stores tRNS chunk if present, and is used to check chunk
        # ordering constraints.
        self.trns = None
        # Stores sBIT chunk if present.
        self.sbit = None

  

    def read(self, lenient=False):
        """
        Read the PNG file and decode it.
        Returns (`width`, `height`, `rows`, `info`).

        May use excessive memory.

        `rows` is a sequence of rows;
        each row is a sequence of values.

        If the optional `lenient` argument evaluates to True,
        checksum failures will raise warnings rather than exceptions.
        """

        def iteridat():
            """Iterator that yields all the ``IDAT`` chunks as strings."""
            while True:
                type, data = self.chunk(lenient=lenient)
                if type == b'IEND':
                    # http://www.w3.org/TR/PNG/#11IEND
                    break
                if type != b'IDAT':
                    continue
                # type == b'IDAT'
                # http://www.w3.org/TR/PNG/#11IDAT
                if self.colormap and not self.plte:
                    warnings.warn("PLTE chunk is required before IDAT chunk")
                yield data

        self.preamble(lenient=lenient)
        raw = decompress(iteridat())

        if self.interlace:
            def rows_from_interlace():
                """Yield each row from an interlaced PNG."""
                # It's important that this iterator doesn't read
                # IDAT chunks until it yields the first row.
                bs = bytearray(itertools.chain(*raw))
                arraycode = 'BH'[self.bitdepth > 8]
                # Like :meth:`group` but
                # producing an array.array object for each row.
                values = self._deinterlace(bs)
                vpr = self.width * self.planes
                for i in range(0, len(values), vpr):
                    row = array(arraycode, values[i:i+vpr])
                    yield row
            rows = rows_from_interlace()
        else:
            rows = self._iter_bytes_to_values(self._iter_straight_packed(raw))
        info = dict()
        for attr in 'greyscale alpha planes bitdepth interlace'.split():
            info[attr] = getattr(self, attr)
        info['size'] = (self.width, self.height)
        for attr in 'gamma transparent background'.split():
            a = getattr(self, attr, None)
            if a is not None:
                info[attr] = a
        if getattr(self, 'x_pixels_per_unit', None):
            info['physical'] = Resolution(self.x_pixels_per_unit,
                                          self.y_pixels_per_unit,
                                          self.unit_is_meter)
        if self.plte:
            info['palette'] = self.palette()
        return self.width, self.height, rows, info

    def read_flat(self):
        """
        Read a PNG file and decode it into a single array of values.
        Returns (*width*, *height*, *values*, *info*).

        May use excessive memory.

        `values` is a single array.

        The :meth:`read` method is more stream-friendly than this,
        because it returns a sequence of rows.
        """

        x, y, pixel, info = self.read()
        arraycode = 'BH'[info['bitdepth'] > 8]
        pixel = array(arraycode, itertools.chain(*pixel))
        self.file.close() #je clos le fichier avant de retourner les valeurs
        return x, y, pixel, info


def decompress(data_blocks):
    """
    `data_blocks` should be an iterable that
    yields the compressed data (from the ``IDAT`` chunks).
    This yields decompressed byte strings.
    """

    # Currently, with no max_length parameter to decompress,
    # this routine will do one yield per IDAT chunk: Not very
    # incremental.
    d = zlib.decompressobj()
    # Each IDAT chunk is passed to the decompressor, then any
    # remaining state is decompressed out.
    for data in data_blocks:
        # :todo: add a max_length argument here to limit output size.
        yield bytearray(d.decompress(data))
    yield bytearray(d.flush())


def check_bitdepth_colortype(bitdepth, colortype):
    """
    Check that `bitdepth` and `colortype` are both valid,
    and specified in a valid combination.
    Returns (None) if valid, raise an Exception if not valid.
    """

    if bitdepth not in (1, 2, 4, 8, 16):
        raise FormatError("invalid bit depth %d" % bitdepth)
    if colortype not in (0, 2, 3, 4, 6):
        raise FormatError("invalid colour type %d" % colortype)
    # Check indexed (palettized) images have 8 or fewer bits
    # per pixel; check only indexed or greyscale images have
    # fewer than 8 bits per pixel.
    if colortype & 1 and bitdepth > 8:
        raise FormatError(
            "Indexed images (colour type %d) cannot"
            " have bitdepth > 8 (bit depth %d)."
            " See http://www.w3.org/TR/2003/REC-PNG-20031110/#table111 ."
            % (bitdepth, colortype))
    if bitdepth < 8 and colortype not in (0, 3):
        raise FormatError(
            "Illegal combination of bit depth (%d)"
            " and colour type (%d)."
            " See http://www.w3.org/TR/2003/REC-PNG-20031110/#table111 ."
            % (bitdepth, colortype))


def undo_filter_sub(filter_unit, scanline, previous, result):
    """Undo sub filter."""

    ai = 0
    # Loops starts at index fu.  Observe that the initial part
    # of the result is already filled in correctly with
    # scanline.
    for i in range(filter_unit, len(result)):
        x = scanline[i]
        a = result[ai]
        result[i] = (x + a) & 0xff
        ai += 1


def undo_filter_up(filter_unit, scanline, previous, result):
    """Undo up filter."""

    for i in range(len(result)):
        x = scanline[i]
        b = previous[i]
        result[i] = (x + b) & 0xff


def undo_filter_average(filter_unit, scanline, previous, result):
    """Undo up filter."""

    ai = -filter_unit
    for i in range(len(result)):
        x = scanline[i]
        if ai < 0:
            a = 0
        else:
            a = result[ai]
        b = previous[i]
        result[i] = (x + ((a + b) >> 1)) & 0xff
        ai += 1


def undo_filter_paeth(filter_unit, scanline, previous, result):
    """Undo Paeth filter."""

    # Also used for ci.
    ai = -filter_unit
    for i in range(len(result)):
        x = scanline[i]
        if ai < 0:
            a = c = 0
        else:
            a = result[ai]
            c = previous[ai]
        b = previous[i]
        p = a + b - c
        pa = abs(p - a)
        pb = abs(p - b)
        pc = abs(p - c)
        if pa <= pb and pa <= pc:
            pr = a
        elif pb <= pc:
            pr = b
        else:
            pr = c
        result[i] = (x + pr) & 0xff
        ai += 1

#####fin de la partie dedier a la lecture du png


class GcodeExport(inkex.Effect):

  def __init__(self):
    """init the effetc library and get options from gui"""
    inkex.Effect.__init__(self)
    
    # Options du menu
    self.arg_parser.add_argument("-d", "--directory", dest="directory", default="/home/",help="Directory for files") ####check_dir
    self.arg_parser.add_argument("-f", "--filename",  dest="filename", default="-1.0", help="File name")
    self.arg_parser.add_argument("--add-numeric-suffix-to-filename",  type=inkex.Boolean, dest="add_numeric_suffix_to_filename", default=True,help="Add numeric suffix to filename")
    self.arg_parser.add_argument("--bg_color",dest="bg_color",default="",help="")
    self.arg_parser.add_argument("--resolution", type=int, dest="resolution", default="5",help="") #Usare il valore su float(xy)/resolution e un case per i DPI dell export
    # passage en gris
    self.arg_parser.add_argument("--grayscale_type", type=int, dest="grayscale_type", default="1",help="")
    # Noir et blanc
    self.arg_parser.add_argument("--conversion_type", type=int, dest="conversion_type", default="1",help="")
    # Options lier au noir et blanc
    self.arg_parser.add_argument("--BW_threshold", type=int, dest="BW_threshold", default="128",help="")
    self.arg_parser.add_argument("--grayscale_resolution", type=int, dest="grayscale_resolution", default="1",help="")
    #Vitesse de déplacement
    self.arg_parser.add_argument("--speed_ON", type=int, dest="speed_ON", default="200",help="")
    # Miroir sur l'axe Y
    self.arg_parser.add_argument("--flip_y", type=inkex.Boolean, dest="flip_y", default=False,help="")
    # cycle du home ou pas
    self.arg_parser.add_argument("--homing", type=int, dest="homing", default="1",help="")
    # Commands
    self.arg_parser.add_argument("--Power_ON", type=int, dest="Power_ON", default="3000",help="")
    self.arg_parser.add_argument("--laseron",  dest="laseron", default="M4", help="M4 for GRBL M104 for Marlin")
    self.arg_parser.add_argument("--laseroff", dest="laseroff", default="M5", help="M5 for GRBL M105 for Marlin")
    self.arg_parser.add_argument("--laserpause",  dest="laserpause", default="G4 P0", help="")
    
    # Prévisualisation, ouvre une fenêtre avec juste la conversion de l'image
    self.arg_parser.add_argument("--preview_only", type=inkex.Boolean, dest="preview_only", default=False,help="")

  def effect(self):
    
    current_file = self.options.input_file
    bg_color = self.options.bg_color
     
    #Check directory
    
    if (os.path.isdir(self.options.directory)) == True:					
      
      #inkex.errormsg("OK") #DEBUG

      #je rajoute une extension numérique si besoin
      if self.options.add_numeric_suffix_to_filename :
        dir_list = os.listdir(self.options.directory) #je charge les fichierdu répertoire
        temp_name =  self.options.filename
        max_n = 0
        for s in dir_list :
          r = re.match(r"^%s_0*(\d+)%s$"%(re.escape(temp_name),'.png' ), s)
          if r :
            max_n = max(max_n,int(r.group(1)))	
        self.options.filename = temp_name + "_" + ( "0"*(4-len(str(max_n+1))) + str(max_n+1) ) #je reconstruit le nom s il existe déja avec les options
      
      suffix = ""
      if self.options.conversion_type == 1:
        suffix = "_BWfix_"+str(self.options.BW_threshold)+"_"
      elif self.options.conversion_type == 2:
        suffix = "_BWrnd_"
      else:
        if self.options.grayscale_resolution == 1:
          suffix = "_Gray_256_"
        elif self.options.grayscale_resolution == 2:
          suffix = "_Gray_128_"
        elif self.options.grayscale_resolution == 4:
          suffix = "_Gray_64_"
        elif self.options.grayscale_resolution == 8:
          suffix = "_Gray_32_"
        elif self.options.grayscale_resolution == 16:
          suffix = "_Gray_16_"
        elif self.options.grayscale_resolution == 32:
          suffix = "_Gray_8_"
        else:
          suffix = "_Gray_"
        
      
      pos_file_png_exported = os.path.join(self.options.directory,self.options.filename+".png") 
      pos_file_png_BW = os.path.join(self.options.directory,self.options.filename+suffix+"preview.png") 
      pos_file_gcode = os.path.join(self.options.directory,self.options.filename+suffix+".gcode") 
      

      #j'exporte l'image et je lance la convertion
      self.exportPage(pos_file_png_exported,current_file,bg_color)

      self.PNGtoGcode(pos_file_png_exported,pos_file_png_BW,pos_file_gcode)
            
      
    else:
      inkex.errormsg("Directory does not exist! Please specify existing directory!")
            
    
  def exportPage(self,pos_file_png_exported,current_file,bg_color):		
     #l'export utilise simplement la ligne de commande de inskcape
    if self.options.resolution == 1:
      DPI = 25.4
    elif self.options.resolution == 2:
      DPI = 50.8
    elif self.options.resolution == 5:
      DPI = 127
    else:
      DPI = 254
    
    #inkscape --export-type="png" my_file.svg
    #--export-filename
    # old pour la version 0.92 command="inkscape -C -e \"%s\" -b\"%s\" %s -d %s" % (pos_file_png_exported,bg_color,current_file,DPI) #Comando da linea di comando per esportare in PNG
    #inkscape -C -e "/home/David/test_0001.png" -b"#ffffff" /tmp/ink_ext_XXXXXX.svg5H39R1 -d 25.4
    command='inkscape --export-type="png" --export-filename=%s %s -b \"%s\"  -d %s' % (pos_file_png_exported,current_file,bg_color,DPI)
    #inkex.errormsg(command)
    #pour le passage en version python 3 je l'encadre pour être sur que le process ce ferme bien
    with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
        _ = p.wait()

  def PNGtoGcode(self,pos_file_png_exported,pos_file_png_BW,pos_file_gcode):
    #on transforme une image en gcode

    # inkex.errormsg(pos_file_png_exported)
    MyPNG = Reader(pos_file_png_exported) #File PNG generato
    w, h, pixels, metadata = MyPNG.read_flat()
    
    # inkex.errormsg(str(w)+ " " + str(h)+" " + str(metadata))
    # inkex.errormsg(pixels)
   
    matrice = [[255 for i in range(w)]for j in range(h)]  #je récupère l'array
    ######## GENERO IMMAGINE IN BIANCO E NERO ########
    #niveau pour le noir et blanc
    B=255
    N=0 
    
    #et on le traite pixel par pixel
    for y in range(h): # y varia da 0 a h-1
      for x in range(w): # x varia da 0 a w-1
        pixel_position = (x + y * w)*4 if metadata['alpha'] else (x + y * w)*3
        
        #options de passage en niveau de gris
        if self.options.grayscale_type == 1:
            couleur=int(pixels[pixel_position]*0.21 + pixels[(pixel_position+1)]*0.71 + pixels[(pixel_position+2)]*0.07)
        elif self.options.grayscale_type == 2:
            couleur = int((pixels[pixel_position] + pixels[(pixel_position+1)]+ pixels[(pixel_position+2)]) / 3 )	 
        elif self.options.grayscale_type == 3:
            couleur = int(pixels[pixel_position]) #Red
        elif self.options.grayscale_type == 4:
            couleur = int(pixels[(pixel_position+1)]) # green	
        elif self.options.grayscale_type == 5:
            couleur = int(pixels[(pixel_position+2)]) # blue
        elif self.options.grayscale_type == 6:
            couleur= int(max(pixels[pixel_position] , pixels[(pixel_position+1)] , pixels[(pixel_position+2)]))
        elif self.options.grayscale_type == 7:
            couleur=int(min(pixels[pixel_position] , pixels[(pixel_position+1)] , pixels[(pixel_position+2)]))
           
        # et/ou transformation en noir et blanc
        # couleur = int(min(pixels[pixel_position] , pixels[(pixel_position+1)] , pixels[(pixel_position+2)]))
        if self.options.conversion_type == 1:
            if couleur >= self.options.BW_threshold:
              couleur = N
            else:
              couleur = B
        elif self.options.conversion_type == 2:
            if couleur >= randint(self.options.BW_threshold-50,self.options.BW_threshold+50):
              couleur = N
            else:
              couleur = B     
                  
        if couleur <= 1:
          couleur = 0
        if couleur >= 254:
          couleur = 255

        matrice[y][x]=couleur
             
           
    #### GENERO IL FILE GCODE ####
    if self.options.preview_only == False: #Genero Gcode solo se devo
    
      if self.options.flip_y == False: #Inverto asse Y solo se flip_y = False     
        #-> coordinate Cartesiane (False) Coordinate "informatiche" (True)
        matrice.reverse()				


      F_G1 = self.options.speed_ON
      Scala = self.options.resolution
      
      with open(pos_file_gcode, 'w') as file_gcode: #Creo il file
          
          gcode=""
          #Configurazioni iniziali standard Gcode
          gcode+='( Generated with: Raster 2 Laser Gcode generator)\n'

          #bCNC and GRBL compatibility
          gcode+='M5 S0\n'
          gcode+='G90\n'
          gcode+='G21\n'
          gcode+='G0 ' +' F' + str(F_G1) + '\n' #meme vitesse que la gravure pour eviter les accoups
          gcode+='G1 ' +' F' + str(F_G1) + '\n' #meme vitesse que la gravure pour eviter les accoups
          gcode+=self.options.laseron + ' S0\n'
          #partie creation du gcode

          # for y in range(h):
          #    matrice[y].append(B)
          # w = w+1
          
            #noir et blanc
          #if self.options.conversion_type != 3: 
          for y in range(h):
          # if y+1 <h:
          #     if matrice[y][0] != matrice[y+1][0]:
              gcode+='G0 Y'+ str(float(y)/Scala) + ' S0\n'
              for x in range(w):
                  check_pixel=False
                  if y % 2 == 0 : #je check pour savoir si le pixel suivant est identique
                      dir_x=x
                      if dir_x > 1:
                          check_pixel=(matrice[y][dir_x-1] != matrice[y][dir_x])
                  else:
                      dir_x=(w-1)-x
                      if dir_x < (w-1):
                          check_pixel=(matrice[y][dir_x+1] != matrice[y][dir_x])
        
                  #je ne travaille que sur des valeur non blanche
                  if matrice[y][dir_x] !=N:
                      Type_G='G0 X'
                  else:
                      Type_G='G1 X'
                      
                  if check_pixel:#je ne rajoute une ligne que si j'ai une modification de la valeur de couleur
                      gcode+=Type_G + str(float(dir_x)/Scala) + ' S'+ str(round((255 - matrice[y][dir_x])/255*self.options.Power_ON)) +'\n'
 

          #Configurazioni finali standard Gcode
          gcode+=self.options.laseroff  +'\n'
          gcode+='G0 X0 Y0; home\n'
          #HOMING
          if self.options.homing == 1:
            gcode+='G28\n'
          elif self.options.homing == 2:
            gcode+='$H\n'
          else:
            pass
          
          file_gcode.write(gcode)
          #file_gcode.close() #Chiudo il file
      

if __name__=="__main__":
   GcodeExport().run()
   
   exit()


