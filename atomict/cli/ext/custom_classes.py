import click


class DefaultCommandGroup(click.Group):
    """Allow a default command for a group.
    
    This class extends Click's Group to support a default command that gets executed
    when no subcommand is specified or when arguments don't match any known command.
    """

    def command(self, *args, **kwargs):
        default_command = kwargs.pop('default_command', False)
        if default_command and not args:
            kwargs['name'] = kwargs.get('name', ' ')
        decorator = super(DefaultCommandGroup, self).command(*args, **kwargs)

        if default_command:
            def new_decorator(f):
                cmd = decorator(f)
                self.default_command = cmd.name
                return cmd

            return new_decorator

        return decorator

    def resolve_command(self, ctx, args):
        try:
            # test if the command parses
            return super(DefaultCommandGroup, self).resolve_command(ctx, args)
        except click.UsageError:
            # command did not parse, assume it is the default command
            args.insert(0, self.default_command)
            return super(DefaultCommandGroup, self).resolve_command(ctx, args)
