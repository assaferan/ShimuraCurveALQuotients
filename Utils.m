// General utility intrinsics.

intrinsic WriteStderr(s::MonStgElt)
{ write to stderr }
  E := Open("/dev/stderr", "a");
  Write(E, s);
  Flush(E);
end intrinsic;

intrinsic WriteStderr(e::Err)
{ write to stderr }
  WriteStderr(Sprint(e) cat "\n");
end intrinsic;
